#!/bin/bash
# Publish the FUSE container image to the GitHub Container Registry as a
# multi-arch tag:
#
#     ghcr.io/projecttorreypines/fuse:<version>        (manifest list)
#         ├── ghcr.io/projecttorreypines/fuse:<version>-amd64
#         └── ghcr.io/projecttorreypines/fuse:<version>-arm64
#
# The version tag (and `latest`) is a manifest list: one tag that points at
# both per-arch images, so `docker run ghcr.io/projecttorreypines/fuse:latest`
# automatically pulls the image matching the client's architecture.
#
# Publishing is therefore two steps:
#
#   1. Per-arch push — run AFTER build.sh, on the same node (the image lives
#      in the node-local podman store). Pushes to :<version>-<ARCH> and does
#      NOT touch :<version> or :latest.
#
#          FUSE_ENVIRONMENT=v1.2.0 ./deploy/omega-container/publish_ghcr.sh
#
#      ARCH defaults to amd64 (the omega build). The arm64 image is built and
#      pushed from an Apple Silicon machine with docker (see the README).
#
#   2. Manifest publish — run from ANY machine once BOTH per-arch tags exist.
#      Combines them into :<version> and :latest. Only do this after both
#      arch images passed their acceptance tests: this is the step that moves
#      `latest` for all users.
#
#          MANIFEST=1 FUSE_ENVIRONMENT=v1.2.0 ./deploy/omega-container/publish_ghcr.sh
#
# NOTE: the version comes from the FUSE_ENVIRONMENT environment variable (or,
# when unset, the latest GitHub release) — NOT from a positional argument.
#
# Requires the `gh` CLI authenticated with the write:packages scope:
#
#     gh auth refresh --hostname github.com -s write:packages
#
# Notes:
# - The amd64 image built with the default (universal) JULIA_CPU_TARGET runs
#   on NERSC Perlmutter (znver3), GA omega (cascadelake/znver2/znver3), and
#   any x86_64 laptop (generic fallback). The arm64 image ships without a
#   sysimage (no aarch64 FUSE sysimage can link) and serves Apple Silicon.
# - The FIRST push creates the package with PRIVATE visibility. To let other
#   sites pull without authentication, make it public once (org owners /
#   package admins) at:
#   https://github.com/orgs/ProjectTorreyPines/packages/container/package/fuse
#   -> Package settings -> Change visibility.

set -euo pipefail

if [[ $# -gt 0 ]]; then
    echo "ERROR: this script takes no positional arguments. Set the version" >&2
    echo "       via FUSE_ENVIRONMENT=$1 instead." >&2
    exit 1
fi

if [[ -n "${FUSE_ENVIRONMENT:-}" ]]; then
    version="$FUSE_ENVIRONMENT"
else
    version="$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)"
fi
if [[ -z "$version" || "$version" == "null" ]]; then
    echo "ERROR: could not determine FUSE version (set FUSE_ENVIRONMENT)." >&2
    exit 1
fi

dest="ghcr.io/projecttorreypines/fuse"

if ! command -v gh >/dev/null 2>&1 || ! gh auth status >/dev/null 2>&1; then
    echo "ERROR: gh CLI not authenticated. Run 'gh auth login' and add the" >&2
    echo "       write:packages scope (see header)." >&2
    exit 1
fi
user="$(gh api /user --jq .login)"

# Node-local podman storage (omega worker nodes): the default graphRoot is in
# NFS $HOME, which is slow and eats quota. Fall back to podman defaults on
# machines without /local-scratch (e.g. a laptop running the MANIFEST step).
podman=(podman)
scratch="/local-scratch/$USER"
if [[ -d "$scratch" ]]; then
    podman=(podman --root "$scratch/podman-storage" --runroot "$scratch/podman-run")
fi

if [[ "${MANIFEST:-0}" == "1" ]]; then
    # --- step 2: combine the per-arch tags into :<version> and :latest ------
    if command -v docker >/dev/null 2>&1 && docker buildx version >/dev/null 2>&1; then
        # Registry-side manifest creation — copies no blobs.
        echo "### Logging in to ghcr.io as $user (docker)"
        gh auth token | docker login ghcr.io -u "$user" --password-stdin
        trap 'docker logout ghcr.io >/dev/null 2>&1 || true' EXIT
        echo "### Creating manifest $dest:$version (+ latest) from -amd64 + -arm64"
        docker buildx imagetools create \
            -t "$dest:$version" -t "$dest:latest" \
            "$dest:$version-amd64" "$dest:$version-arm64"
        docker buildx imagetools inspect "$dest:$version"
    else
        echo "### Logging in to ghcr.io as $user (podman)"
        gh auth token | "${podman[@]}" login ghcr.io -u "$user" --password-stdin
        trap '"${podman[@]}" manifest rm "fuse-manifest-$version" >/dev/null 2>&1 || true
              "${podman[@]}" logout ghcr.io >/dev/null 2>&1 || true' EXIT
        echo "### Creating manifest $dest:$version (+ latest) from -amd64 + -arm64"
        "${podman[@]}" manifest rm "fuse-manifest-$version" >/dev/null 2>&1 || true
        "${podman[@]}" manifest create "fuse-manifest-$version"
        # --all is REQUIRED: without it, `manifest add` tries to select the
        # image matching the *host* architecture from the remote reference and
        # fails on the foreign-arch tag ("no image found in image index for
        # architecture amd64"). --all copies every entry, which also carries
        # over the buildx attestation manifest. (verified on podman 4.9.4)
        "${podman[@]}" manifest add --all "fuse-manifest-$version" "docker://$dest:$version-amd64"
        "${podman[@]}" manifest add --all "fuse-manifest-$version" "docker://$dest:$version-arm64"
        "${podman[@]}" manifest push --all "fuse-manifest-$version" "docker://$dest:$version"
        "${podman[@]}" manifest push --all "fuse-manifest-$version" "docker://$dest:latest"
        # No docker:// prefix here: `manifest add` accepts that transport but
        # `manifest inspect` does not ("unsupported transport docker for looking
        # up local images") — it takes a plain remote reference.
        "${podman[@]}" manifest inspect "$dest:$version"
    fi
    echo
    echo "### Published $dest:$version and $dest:latest (multi-arch)"
    exit 0
fi

# --- step 1: push the locally-built image to its per-arch tag ---------------
arch="${ARCH:-amd64}"
image="localhost/fuse:$version"

if ! "${podman[@]}" image exists "$image"; then
    echo "ERROR: $image not found in the podman store on this node." >&2
    echo "       Run build.sh here first (the store is node-local)." >&2
    exit 1
fi

echo "### Logging in to ghcr.io as $user"
gh auth token | "${podman[@]}" login ghcr.io -u "$user" --password-stdin

# Always log out afterwards so no registry credential lingers on the node.
trap '"${podman[@]}" logout ghcr.io >/dev/null 2>&1 || true' EXIT

echo "### Pushing $image -> $dest:$version-$arch (several GB, takes a while)"
"${podman[@]}" push "$image" "$dest:$version-$arch"

echo
echo "### Published $dest:$version-$arch"
echo "Users pull :$version / :latest, which are manifest lists — publish those"
echo "with MANIFEST=1 once BOTH -amd64 and -arm64 tags exist:"
echo "    MANIFEST=1 FUSE_ENVIRONMENT=$version $0"
