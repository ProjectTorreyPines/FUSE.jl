#!/bin/bash
# Publish the built FUSE container image to the GitHub Container Registry:
#
#     ghcr.io/projecttorreypines/fuse:<version>
#
# Run AFTER build.sh, on the same node (the image lives in the node-local
# podman store). Requires the `gh` CLI authenticated with the write:packages
# scope:
#
#     gh auth refresh --hostname github.com -s write:packages
#
# Usage:
#     FUSE_ENVIRONMENT=v1.1.5 ./deploy/omega-container/publish_ghcr.sh
#
# Optional:
#     LATEST=1   also tag/push ghcr.io/projecttorreypines/fuse:latest
#
# Notes:
# - The image built with the default (universal) JULIA_CPU_TARGET runs on
#   NERSC Perlmutter (znver3), GA omega (cascadelake/znver2/znver3), and any
#   x86_64 laptop (generic fallback), so one published tag serves all sites.
# - The FIRST push creates the package with PRIVATE visibility. To let other
#   sites pull without authentication, make it public once (org owners /
#   package admins) at:
#   https://github.com/orgs/ProjectTorreyPines/packages/container/package/fuse
#   -> Package settings -> Change visibility.

set -euo pipefail

if [[ -n "${FUSE_ENVIRONMENT:-}" ]]; then
    version="$FUSE_ENVIRONMENT"
else
    version="$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)"
fi
if [[ -z "$version" || "$version" == "null" ]]; then
    echo "ERROR: could not determine FUSE version (set FUSE_ENVIRONMENT)." >&2
    exit 1
fi

scratch="/local-scratch/$USER"
podman=(podman --root "$scratch/podman-storage" --runroot "$scratch/podman-run")
image="localhost/fuse:$version"
dest="ghcr.io/projecttorreypines/fuse"

if ! "${podman[@]}" image exists "$image"; then
    echo "ERROR: $image not found in the podman store on this node." >&2
    echo "       Run build.sh here first (the store is node-local)." >&2
    exit 1
fi
if ! command -v gh >/dev/null 2>&1 || ! gh auth status >/dev/null 2>&1; then
    echo "ERROR: gh CLI not authenticated. Run 'gh auth login' and add the" >&2
    echo "       write:packages scope (see header)." >&2
    exit 1
fi

user="$(gh api /user --jq .login)"
echo "### Logging in to ghcr.io as $user"
gh auth token | "${podman[@]}" login ghcr.io -u "$user" --password-stdin

# Always log out afterwards so no registry credential lingers on the node.
trap '"${podman[@]}" logout ghcr.io >/dev/null 2>&1 || true' EXIT

echo "### Pushing $image -> $dest:$version (several GB, takes a while)"
"${podman[@]}" push "$image" "$dest:$version"

if [[ "${LATEST:-0}" == "1" ]]; then
    echo "### Pushing $dest:latest"
    "${podman[@]}" push "$image" "$dest:latest"
fi

echo
echo "### Published $dest:$version"
echo "If this was the first push, the package is PRIVATE — make it public in"
echo "the package settings so all sites can pull it unauthenticated:"
echo "    https://github.com/orgs/ProjectTorreyPines/packages/container/package/fuse"
