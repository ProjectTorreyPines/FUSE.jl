#!/bin/bash
# One-command FUSE container setup for NERSC Perlmutter Jupyter.
#
# Gets the published FUSE image usable on this account (reusing the m3739
# shared squash store when possible, otherwise pulling from ghcr.io) and
# installs the Jupyter kernelspec. FUSE then runs from
# https://jupyter.nersc.gov — no Julia install, no compilation.
#
#   bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/deploy/perlmutter-container/install_fuse_container_nersc.sh)
#
# Optional:
#   FUSE_ENVIRONMENT=v1.2.0   image version (default: latest FUSE.jl release)
#   THREADS=8                 Julia threads for the kernel (default 1)
#   SQUASH_DIR=...            force a specific shared image store
#
# This is the NERSC counterpart of omega's
# `module load fuse-container && fuse-container lab`: JupyterHub plays the role
# of the host-side JupyterLab, so setup ends in the browser rather than by
# spawning a server here.

set -euo pipefail

shared_store="/global/cfs/cdirs/m3739/shared_images"
registry_image_base="ghcr.io/projecttorreypines/fuse"

command -v podman-hpc >/dev/null || { echo "ERROR: podman-hpc not found — run this on Perlmutter." >&2; exit 1; }

if [[ -n "${FUSE_ENVIRONMENT:-}" ]]; then
    version="$FUSE_ENVIRONMENT"
else
    version="$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)"
fi
[[ -n "$version" && "$version" != "null" ]] || { echo "ERROR: could not determine FUSE version (set FUSE_ENVIRONMENT)." >&2; exit 1; }

has_image() {  # has_image <extra podman-hpc args...>
    podman-hpc "$@" images --format '{{.Repository}}:{{.Tag}}' 2>/dev/null \
        | grep -qx "localhost/fuse:$version"
}

# --- pick an image store -----------------------------------------------------
# Prefer the m3739 shared store when it already holds the requested version;
# fall back to pulling from ghcr.io into the per-user store when the shared
# store is stale (a newer release is out) or not accessible.
squash_dir="${SQUASH_DIR:-}"
if [[ -z "$squash_dir" ]]; then
    if ! id -nG | tr ' ' '\n' | grep -qx m3739 || [[ ! -r "$shared_store" ]]; then
        echo "### m3739 shared store not accessible — using your per-user image store"
    elif has_image --squash-dir "$shared_store"; then
        squash_dir="$shared_store"
        echo "### Using the shared m3739 image store (no pull needed): $squash_dir"
    else
        echo "### Shared m3739 store does not have fuse:$version yet — falling back to a registry pull"
    fi
fi

if [[ -n "$squash_dir" ]]; then
    has_image --squash-dir "$squash_dir" \
        || { echo "ERROR: localhost/fuse:$version not found in $squash_dir" >&2; exit 1; }
elif ! has_image; then
    echo "### Pulling $registry_image_base:$version (~20 GB, several minutes)"
    podman-hpc pull "$registry_image_base:$version"
    podman-hpc tag "$registry_image_base:$version" "localhost/fuse:$version"
    podman-hpc migrate "fuse:$version"
else
    echo "### localhost/fuse:$version already in your image store"
fi

# --- install the Jupyter kernelspec ------------------------------------------
# Reuse the canonical generator. When this script runs from a FUSE checkout the
# generator sits next to it; when run via `bash <(curl ...)` there is no
# checkout, so fetch the generator and its template from the same branch.
scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || echo /dev/fd)"
if [[ -f "$scriptdir/install_kernel.sh" && -f "$scriptdir/kernel.json.template" ]]; then
    installer="$scriptdir/install_kernel.sh"
else
    raw="https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/deploy/perlmutter-container"
    tmp="$(mktemp -d)"
    trap 'rm -rf "$tmp"' EXIT
    curl -fsSL "$raw/install_kernel.sh" -o "$tmp/install_kernel.sh"
    curl -fsSL "$raw/kernel.json.template" -o "$tmp/kernel.json.template"
    chmod +x "$tmp/install_kernel.sh"
    installer="$tmp/install_kernel.sh"
fi

SQUASH_DIR="$squash_dir" FUSE_ENVIRONMENT="$version" THREADS="${THREADS:-1}" "$installer"

echo
echo "### FUSE container ready."
echo "Open https://jupyter.nersc.gov, start a 'Login Node' server, and select"
echo "the 'Julia FUSE-$version' kernel. First cell to try:"
echo "    using FUSE; ini, act = FUSE.case_parameters(:D3D, :L_mode); dd = FUSE.init(ini, act)"
