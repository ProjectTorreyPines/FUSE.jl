#!/bin/bash
# Install a host-side Jupyter kernelspec that runs FUSE inside the Singularity
# container on omega. Jupyter discovers kernels under
# $HOME/.local/share/jupyter/kernels/. The kernel's argv runs the in-container
# Julia (with the FUSE sysimage) through the `fuse-container` launcher, which
# manages the squashfuse mount so it is cleaned up even when Jupyter kills the
# kernel abruptly (see fuse-container header). singularity's default binds
# ($HOME, /tmp) let the kernel reach the connection file.
#
# Run this AFTER ./build.sh has produced the SIF. Usage:
#   module load singularity/3.11.3
#   SIF=/fusion/.../fuse_v2.6.0.sif ./deploy/omega-container/install_kernel.sh
#
# Optional:
#   THREADS=8  number of Julia threads the kernel starts with (default 1)

set -euo pipefail

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${SIF:?set SIF=/path/to/fuse_<version>.sif}"
if [[ ! -r "$SIF" ]]; then
    echo "ERROR: SIF not found/readable: $SIF" >&2
    exit 1
fi
threads="${THREADS:-1}"

# Version for the kernel name: parse from the SIF filename (fuse_<version>.sif).
version="$(basename "$SIF" .sif)"
version="${version#fuse_}"

if ! command -v singularity >/dev/null 2>&1; then
    echo "ERROR: singularity not found. Run 'module load singularity/3.11.3' first." >&2
    exit 1
fi
# Kernel argv cannot 'module load', so bake the absolute singularity bin dir
# onto the kernel's PATH (the launcher needs `singularity`).
singularity_dir="$(dirname "$(command -v singularity)")"

# Prefer the launcher published next to the SIF (stable path), else the repo
# copy. squashfuse is found by the launcher itself from <dir-of-SIF>/bin.
if [[ -x "$(dirname "$SIF")/fuse-container" ]]; then
    launcher="$(dirname "$SIF")/fuse-container"
else
    launcher="$scriptdir/fuse-container"
fi
if [[ ! -x "$(dirname "$SIF")/bin/squashfuse" ]] && ! command -v squashfuse >/dev/null 2>&1; then
    echo "WARNING: squashfuse not found next to the SIF (<dir-of-SIF>/bin/) or on" >&2
    echo "         PATH. The kernel will fail to start until it is available." >&2
fi

kernel_dir="$HOME/.local/share/jupyter/kernels/fuse-$version"
mkdir -p "$kernel_dir"

display="Julia FUSE-$version container ($threads thread(s))"

sed -e "s|__SINGULARITY_DIR__|$singularity_dir|g" \
    -e "s|__LAUNCHER__|$launcher|g" \
    -e "s|__SIF__|$SIF|g" \
    -e "s|__THREADS__|$threads|g" \
    -e "s|__DISPLAY__|$display|g" \
    "$scriptdir/kernel.json.template" > "$kernel_dir/kernel.json"

echo
echo "### Installed kernelspec at $kernel_dir/kernel.json"
echo "    launcher: $launcher"
echo "Select the '$display' kernel in Jupyter (see docs/src/install_omega.md"
echo "for the SSH-tunnel workflow), or test it headless with:"
echo "    python3 $scriptdir/test_kernel_headless.py fuse-$version"
