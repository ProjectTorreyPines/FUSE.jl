#!/bin/bash
# Install a host-side Jupyter kernelspec that runs FUSE inside the podman-hpc
# container. NERSC JupyterHub discovers kernels under
# $HOME/.local/share/jupyter/kernels/. The kernel's argv wraps the in-container
# Julia (with the FUSE sysimage) using `podman-hpc run --jupyter`; the
# --jupyter flag bind-mounts /tmp and $HOME so the kernel can connect and write
# notebooks. `--network host` is required: rootless podman otherwise gives the
# container its own network namespace, so the ZMQ ports the kernel binds are
# unreachable from the Jupyter server on the host and the kernel never connects.
#
# Run this AFTER ./build.sh has built and migrated the image. Usage:
#   FUSE_ENVIRONMENT=v1.2.0 ./deploy/perlmutter-container/install_kernel.sh
#
# Optional:
#   THREADS=8      number of Julia threads the kernel starts with (default 1)
#   SQUASH_DIR=... use an image from a shared squash dir instead of the
#                  per-user store. The generated kernel will pass
#                  `podman-hpc --squash-dir <dir> run ...`. Example (project
#                  image shared with all of m3739):
#                    SQUASH_DIR=/global/cfs/cdirs/m3739/shared_images \
#                    FUSE_ENVIRONMENT=v1.2.0 ./install_kernel.sh

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

# Confirm the image is usable before installing the kernelspec. The kernel
# launches IJulia with `-e "import IJulia; IJulia.run_kernel()"` (the modern
# IJulia entrypoint — running src/kernel.jl as a script no longer works).
echo "### Checking $image${squash_dir:+ (squash-dir: $squash_dir)}"
if ! podman-hpc "${podman_global[@]}" run --rm "$image" julia -e 'import IJulia' >/dev/null; then
    echo "ERROR: could not import IJulia from $image." >&2
    echo "       Did you run build.sh (build + migrate) first?" >&2
    exit 1
fi

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

# --jupyter only mounts /tmp and $HOME; also mount the user's scratch (same
# path inside and out, so notebook paths resolve identically) or notebooks and
# data under $PSCRATCH are invisible to the kernel.
if [[ -n "${PSCRATCH:-}" ]]; then
    volume_repl="    \"--volume\",\\n    \"$PSCRATCH:$PSCRATCH\","
else
    volume_repl=""
fi

sed -e "s|__IMAGE__|$image|g" \
    -e "s|__THREADS__|$threads|g" \
    -e "s|__DISPLAY__|$display|g" \
    -e "s|^__SQUASH_ARGS__\$|$squash_repl|" \
    -e "s|^__VOLUME_ARGS__\$|$volume_repl|" \
    "$scriptdir/kernel.json.template" \
  | sed '/^[[:space:]]*$/d' > "$kernel_dir/kernel.json"

# Copy the Julia logos out of the image so JupyterHub shows the Julia tile
# instead of a generic placeholder (IJulia.installkernel does this for
# module-based kernels).
podman-hpc "${podman_global[@]}" run --rm --volume "$kernel_dir:/kout" "$image" \
    julia -e 'import IJulia
              deps = joinpath(dirname(dirname(pathof(IJulia))), "deps")
              for f in readdir(deps)
                  startswith(f, "logo") && cp(joinpath(deps, f), joinpath("/kout", f); force=true)
              end' \
    || echo "WARNING: could not copy kernel logos (kernel still works)"

echo
echo "### Installed kernelspec at $kernel_dir/kernel.json"
[[ -n "$squash_dir" ]] && echo "    (kernel runs the image from squash-dir: $squash_dir)"
echo "Open NERSC JupyterHub (Perlmutter login node) and select the"
echo "'$display' kernel."
