#!/bin/bash
# Install a host-side Jupyter kernelspec that runs FUSE inside the podman-hpc
# container. NERSC JupyterHub discovers kernels under
# $HOME/.local/share/jupyter/kernels/. The kernel's argv wraps the in-container
# Julia (with the FUSE sysimage) using `podman-hpc run --jupyter`; the
# --jupyter flag bind-mounts /tmp and $HOME so the kernel can connect and write
# notebooks.
#
# Run this AFTER ./build.sh has built and migrated the image. Usage:
#   FUSE_ENVIRONMENT=v1.1.3 ./deploy/perlmutter-container/install_kernel.sh
#
# Optional:
#   THREADS=8   number of Julia threads the kernel starts with (default 1)

set -euo pipefail

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -n "${FUSE_ENVIRONMENT:-}" ]]; then
    version="$FUSE_ENVIRONMENT"
else
    version="$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)"
fi
if [[ -z "$version" || "$version" == "null" ]]; then
    echo "ERROR: could not determine FUSE version (set FUSE_ENVIRONMENT)." >&2
    exit 1
fi

threads="${THREADS:-1}"
image="localhost/fuse:$version"

if ! command -v podman-hpc >/dev/null 2>&1; then
    echo "ERROR: podman-hpc not found. Run this on Perlmutter." >&2
    exit 1
fi

# Pull the IJulia kernel.jl path that install_fuse_container.jl recorded in the
# image, so the kernelspec points at the right file inside the container.
echo "### Resolving IJulia kernel path from $image"
kernel_jl="$(podman-hpc run --rm "$image" cat /opt/fuse/ijulia_kernel_path.txt | tr -d '\r\n')"
if [[ -z "$kernel_jl" ]]; then
    echo "ERROR: could not read /opt/fuse/ijulia_kernel_path.txt from $image." >&2
    echo "       Did you run build.sh (build + migrate) first?" >&2
    exit 1
fi
echo "    $kernel_jl"

kernel_dir="$HOME/.local/share/jupyter/kernels/fuse-$version"
mkdir -p "$kernel_dir"

display="Julia FUSE-$version ($threads thread(s))"

sed -e "s|__IMAGE__|$image|g" \
    -e "s|__KERNEL_JL__|$kernel_jl|g" \
    -e "s|__THREADS__|$threads|g" \
    -e "s|__DISPLAY__|$display|g" \
    "$scriptdir/kernel.json.template" > "$kernel_dir/kernel.json"

echo
echo "### Installed kernelspec at $kernel_dir/kernel.json"
echo "Open NERSC JupyterHub (Perlmutter login node) and select the"
echo "'$display' kernel."
