#!/bin/bash
# Build and migrate the FUSE podman-hpc image on Perlmutter.
#
# The sysimage build is heavy (runs the FUSE test suite + tutorial), so run this
# from inside an interactive CPU allocation rather than bare on a login node:
#
#   salloc -N 1 -C cpu -q interactive -A m3739 -t 04:00:00
#   ./deploy/perlmutter-container/build.sh
#
# Override the image tag (defaults to the latest FUSE.jl release) with:
#   FUSE_ENVIRONMENT=v1.1.3 ./deploy/perlmutter-container/build.sh
#
# Share the migrated image across the project (instead of personal SCRATCH):
#   SQUASH_DIR=/global/cfs/cdirs/m3739/shared_images ./build.sh

set -euo pipefail

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# FUSE repo root == build context (Containerfile copies the Makefile from here)
repo_root="$(cd "$scriptdir/../.." && pwd)"

# Resolve the image tag: explicit override, else latest GitHub release name.
if [[ -n "${FUSE_ENVIRONMENT:-}" ]]; then
    version="$FUSE_ENVIRONMENT"
else
    version="$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)"
fi
if [[ -z "$version" || "$version" == "null" ]]; then
    echo "ERROR: could not determine FUSE version (set FUSE_ENVIRONMENT)." >&2
    exit 1
fi

image="fuse:$version"

if ! command -v podman-hpc >/dev/null 2>&1; then
    echo "ERROR: podman-hpc not found. Run this on Perlmutter." >&2
    exit 1
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "WARNING: not inside a Slurm allocation. The sysimage build is heavy;" >&2
    echo "         prefer 'salloc -N 1 -C cpu -q interactive -A m3739 -t 04:00:00' first." >&2
fi

echo "### Building $image (context: $repo_root)"
podman-hpc build -t "$image" -f "$scriptdir/Containerfile" "$repo_root"

echo "### Migrating $image to a squashed read-only image"
if [[ -n "${SQUASH_DIR:-}" ]]; then
    podman-hpc --squash-dir "$SQUASH_DIR" migrate "$image"
else
    podman-hpc migrate "$image"
fi

echo
echo "### Done. Image '$image' is built and migrated."
echo "Run the FUSE REPL with:"
echo "    podman-hpc run --rm -it $image"
echo "Install the Jupyter kernel with:"
echo "    FUSE_ENVIRONMENT=$version ./deploy/perlmutter-container/install_kernel.sh"
