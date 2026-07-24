#!/bin/bash
# Build the FUSE container for GA's omega cluster and export it as a
# Singularity SIF that runs on any omega node.
#
# Pipeline: rootless podman build (node-local storage) -> podman save ->
# singularity build. The SIF is omega's analog of Perlmutter's
# `podman-hpc migrate`: a single read-only file that can live on /fusion
# and be executed anywhere via `singularity run/exec`.
#
# The sysimage build is heavy (runs the FUSE test suite + tutorial, several
# hours), so prefer an interactive allocation. NOTE: omega's `short` partition
# is capped at 30 min — use medium/long:
#
#   salloc -p medium -t 8:00:00 -c 16
#   module load singularity/3.11.3
#   ./deploy/omega-container/build.sh
#
# Override the image tag (defaults to the latest FUSE.jl release):
#   FUSE_ENVIRONMENT=v2.6.0 ./deploy/omega-container/build.sh
#
# Write the SIF to shared storage (instead of node-local /local-scratch),
# so it is visible from every omega node:
#   SIF_DIR=/fusion/ga/projects/ird/ptp/$USER/fuse_containers ./build.sh
#
# Optional:
#   PODMAN_CLEANUP=1  wipe the node-local podman store afterwards (default:
#                     keep it, so rebuilds reuse cached layers)

set -euo pipefail

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# FUSE repo root == build context (Containerfile copies the Makefile from here)
repo_root="$(cd "$scriptdir/../.." && pwd)"

# CPU targets baked into the sysimage. Default: the Containerfile's universal
# set (generic + cascadelake + znver2 + znver3), optimal on omega, Perlmutter,
# and laptops alike — the same image can be published to GHCR for all sites.
# Override with FUSE_CPU_TARGET for a leaner site-specific build.
cpu_target="${FUSE_CPU_TARGET:-generic;cascadelake,-xsaveopt,-rdrnd,clone_all;znver2,-xsaveopt,-rdrnd,clone_all;znver3,-xsaveopt,-rdrnd,clone_all}"

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

if ! command -v podman >/dev/null 2>&1; then
    echo "ERROR: podman not found. Run this on an omega worker node." >&2
    exit 1
fi
if ! command -v singularity >/dev/null 2>&1; then
    echo "ERROR: singularity not found. Run 'module load singularity/3.11.3' first." >&2
    exit 1
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "WARNING: not inside a Slurm allocation. The sysimage build is heavy;" >&2
    echo "         prefer 'salloc -p medium -t 8:00:00 -c 16' first ('short' is" >&2
    echo "         capped at 30 min). podman storage is node-local: run every" >&2
    echo "         step of this script on the same node." >&2
fi

# All container storage goes to node-local /local-scratch: the default podman
# graphRoot and singularity cache sit in NFS $HOME, which is slow for image
# overlays and eats quota. Per-invocation flags — no persistent config files.
scratch="/local-scratch/$USER"
mkdir -p "$scratch/podman-storage" "$scratch/podman-run" \
         "$scratch/singularity-tmp" "$scratch/singularity-cache"
podman=(podman --root "$scratch/podman-storage" --runroot "$scratch/podman-run")
export SINGULARITY_TMPDIR="$scratch/singularity-tmp"
export SINGULARITY_CACHEDIR="$scratch/singularity-cache"

sif_dir="${SIF_DIR:-$scratch}"
mkdir -p "$sif_dir"
sif="$sif_dir/fuse_${version}.sif"

echo "### Building $image (context: $repo_root, CPU target: $cpu_target)"
"${podman[@]}" build -t "$image" \
    --build-arg JULIA_CPU_TARGET="$cpu_target" \
    -f "$scriptdir/../perlmutter-container/Containerfile" "$repo_root"

echo "### Exporting $image -> $sif"
tar="$scratch/fuse_${version}.tar"
"${podman[@]}" save --format docker-archive -o "$tar" "localhost/$image"
singularity build --force "$sif" "docker-archive://$tar"
rm -f "$tar"
chmod a+r "$sif"

if [[ "${PODMAN_CLEANUP:-0}" == "1" ]]; then
    echo "### Cleaning node-local podman store"
    "${podman[@]}" system reset --force
fi

echo
echo "### Done: $sif"
if [[ "$sif_dir" == "$scratch" ]]; then
    echo "NOTE: the SIF is on node-local /local-scratch — copy it to /fusion (or"
    echo "      rebuild with SIF_DIR=...) so other nodes can use it."
fi
echo "Run the FUSE REPL with (put squashfuse on PATH for fast SIF mounting,"
echo "see README.md):"
echo "    singularity run --cleanenv --sif-fuse --bind /fusion $sif"
echo "Install the Jupyter kernel with:"
echo "    SIF=$sif ./deploy/omega-container/install_kernel.sh"
