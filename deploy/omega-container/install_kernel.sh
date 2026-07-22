#!/bin/bash
# Install a host-side Jupyter kernelspec that runs FUSE inside the Singularity
# container on omega. Jupyter discovers kernels under
# $HOME/.local/share/jupyter/kernels/. The kernel's argv wraps the in-container
# Julia (with the FUSE sysimage) in `singularity exec`; singularity's default
# binds ($HOME, /tmp) let the kernel reach the connection file, and
# `--bind /fusion` exposes the shared filesystems.
#
# Run this AFTER ./build.sh has produced the SIF. Usage:
#   module load singularity/3.11.3
#   SIF=/fusion/.../fuse_v2.6.0.sif ./deploy/omega-container/install_kernel.sh
#
# Optional:
#   THREADS=8  number of Julia threads the kernel starts with (default 1)
#
# Fast startup needs squashfuse on PATH so singularity can FUSE-mount the SIF
# instead of extracting it (omega's singularity has no setuid starter). The
# kernel looks for it in <dir-of-SIF>/bin first (the shared containers layout),
# then falls back to whatever is on PATH at install time.

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
# Kernel argv cannot 'module load', so bake the absolute singularity path.
singularity_bin="$(command -v singularity)"

# Locate squashfuse for fast SIF mounting: next to the SIF, or on PATH.
sif_fuse_flag="--sif-fuse"
if [[ -x "$(dirname "$SIF")/bin/squashfuse" ]]; then
    squashfuse_dir="$(dirname "$SIF")/bin"
elif command -v squashfuse >/dev/null 2>&1; then
    squashfuse_dir="$(dirname "$(command -v squashfuse)")"
else
    echo "WARNING: squashfuse not found — the kernel will extract the SIF on" >&2
    echo "         every start (slow). Install squashfuse next to the SIF" >&2
    echo "         (<dir-of-SIF>/bin/squashfuse) and re-run for fast startup." >&2
    squashfuse_dir="/usr/bin"
    sif_fuse_flag=""
fi

kernel_dir="$HOME/.local/share/jupyter/kernels/fuse-$version"
mkdir -p "$kernel_dir"

display="Julia FUSE-$version container ($threads thread(s))"

sed -e "s|__SINGULARITY__|$singularity_bin|g" \
    -e "s|__SQUASHFUSE_DIR__|$squashfuse_dir|g" \
    -e "s|__SIF_FUSE__|$sif_fuse_flag|g" \
    -e "s|__SIF__|$SIF|g" \
    -e "s|__THREADS__|$threads|g" \
    -e "s|__DISPLAY__|$display|g" \
    "$scriptdir/kernel.json.template" > "$kernel_dir/kernel.json"

echo
echo "### Installed kernelspec at $kernel_dir/kernel.json"
echo "Select the '$display' kernel in Jupyter (see docs/src/install_omega.md"
echo "for the SSH-tunnel workflow), or test it headless with:"
echo "    python3 $scriptdir/test_kernel_headless.py fuse-$version"
