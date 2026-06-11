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
#   THREADS=8      number of Julia threads the kernel starts with (default 1)
#   SQUASH_DIR=... use an image from a shared squash dir instead of the
#                  per-user store. The generated kernel will pass
#                  `podman-hpc --squash-dir <dir> run ...`. Example (project
#                  image shared with all of m3739):
#                    SQUASH_DIR=/global/cfs/cdirs/m3739/shared_images \
#                    FUSE_ENVIRONMENT=v1.1.3 ./install_kernel.sh

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
squash_dir="${SQUASH_DIR:-}"

# Global podman-hpc flags applied both when probing the image and in the kernel.
podman_global=()
if [[ -n "$squash_dir" ]]; then
    podman_global=(--squash-dir "$squash_dir")
fi

if ! command -v podman-hpc >/dev/null 2>&1; then
    echo "ERROR: podman-hpc not found. Run this on Perlmutter." >&2
    exit 1
fi

# Pull the IJulia kernel.jl path that install_fuse_container.jl recorded in the
# image, so the kernelspec points at the right file inside the container.
echo "### Resolving IJulia kernel path from $image${squash_dir:+ (squash-dir: $squash_dir)}"
kernel_jl="$(podman-hpc "${podman_global[@]}" run --rm "$image" cat /opt/fuse/ijulia_kernel_path.txt | tr -d '\r\n')"
if [[ -z "$kernel_jl" ]]; then
    echo "ERROR: could not read /opt/fuse/ijulia_kernel_path.txt from $image." >&2
    echo "       Did you run build.sh (build + migrate) first?" >&2
    exit 1
fi
echo "    $kernel_jl"

kernel_dir="$HOME/.local/share/jupyter/kernels/fuse-$version"
mkdir -p "$kernel_dir"

display="Julia FUSE-$version ($threads thread(s))"

# Build the optional "--squash-dir", "<dir>", argv entries (as JSON lines), or
# leave empty so the placeholder line is dropped.
if [[ -n "$squash_dir" ]]; then
    squash_repl="    \"--squash-dir\",\\n    \"$squash_dir\","
else
    squash_repl=""
fi

sed -e "s|__IMAGE__|$image|g" \
    -e "s|__KERNEL_JL__|$kernel_jl|g" \
    -e "s|__THREADS__|$threads|g" \
    -e "s|__DISPLAY__|$display|g" \
    -e "s|^__SQUASH_ARGS__\$|$squash_repl|" \
    "$scriptdir/kernel.json.template" \
  | sed '/^[[:space:]]*$/d' > "$kernel_dir/kernel.json"

echo
echo "### Installed kernelspec at $kernel_dir/kernel.json"
[[ -n "$squash_dir" ]] && echo "    (kernel runs the image from squash-dir: $squash_dir)"
echo "Open NERSC JupyterHub (Perlmutter login node) and select the"
echo "'$display' kernel."
