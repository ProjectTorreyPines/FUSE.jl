#!/usr/bin/env bash
# One-shot FUSE install for laptops (juliaup + Miniconda if needed).
#
# After creating the fuse conda env this activates it, runs
# fusebot install_IJulia (with make / install_ijulia.sh fallbacks), clones
# FuseExamples, and runs the first three cells of fluxmatcher.ipynb.
# Set FUSE_SKIP_VERIFY=1 to skip the notebook solve.
#
# Copy-paste (from any working directory):
#   curl -fsSL https://install.julialang.org | sh -s -- -y && \
#   bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/install_fuse_laptop.sh)
#
# Or, from a FUSE.jl clone:
#   curl -fsSL https://install.julialang.org | sh -s -- -y && bash scripts/install_fuse_laptop.sh

set -euo pipefail

# Companions are fetched from the same scripts/ directory as this file's
# published URL. Override with FUSE_SCRIPT_BASE_URL. Default matches the
# documented one-command install (ProjectTorreyPines/master).
SCRIPT_BASE_URL="${FUSE_SCRIPT_BASE_URL:-https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts}"

# If the caller passed this script's own raw URL as $1 (optional), derive the
# companion base from it so a single curl to any fork/branch stays self-contained:
#   u=https://raw.githubusercontent.com/ORG/FUSE.jl/REF/scripts/install_fuse_laptop.sh
#   bash <(curl -fsSL "$u") "$u"
if [[ "${1:-}" == http://* || "${1:-}" == https://* ]]; then
    SCRIPT_BASE_URL="${FUSE_SCRIPT_BASE_URL:-${1%/*}}"
    shift
fi

resolve_script_dir() {
    local self_dir
    self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -f "${self_dir}/install_fuse_common.sh" ]]; then
        echo "${self_dir}"
        return 0
    fi

    local bundle_dir="${TMPDIR:-/tmp}/fuse-install-$$"
    mkdir -p "${bundle_dir}"
    # Include fluxmatcher verify scripts so the laptop install can finish by
    # running cells 0–2 even when the registry FUSE package is older.
    for file in install_fuse_common.sh install_fuse_julia.jl \
                verify_fluxmatcher_notebook.sh verify_fluxmatcher_notebook.jl; do
        curl -fsSL "${SCRIPT_BASE_URL}/${file}" -o "${bundle_dir}/${file}"
    done
    echo "${bundle_dir}"
}

SCRIPT_DIR="$(resolve_script_dir)"
exec bash "${SCRIPT_DIR}/install_fuse_common.sh" laptop "$@"
