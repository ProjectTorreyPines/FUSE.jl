#!/usr/bin/env bash
# One-shot FUSE install for laptops (juliaup + Miniconda if needed).
#
# Copy-paste (from any working directory):
#   curl -fsSL https://install.julialang.org | sh -s -- -y && \
#   bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/install_fuse_laptop.sh)
#
# Or, from a FUSE.jl clone:
#   curl -fsSL https://install.julialang.org | sh -s -- -y && bash scripts/install_fuse_laptop.sh

set -euo pipefail

SCRIPT_BASE_URL="${FUSE_SCRIPT_BASE_URL:-https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts}"

resolve_script_dir() {
    local self_dir
    self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -f "${self_dir}/install_fuse_common.sh" ]]; then
        echo "${self_dir}"
        return 0
    fi

    local bundle_dir="${TMPDIR:-/tmp}/fuse-install-$$"
    mkdir -p "${bundle_dir}"
    for file in install_fuse_common.sh install_fuse_julia.jl; do
        curl -fsSL "${SCRIPT_BASE_URL}/${file}" -o "${bundle_dir}/${file}"
    done
    echo "${bundle_dir}"
}

SCRIPT_DIR="$(resolve_script_dir)"
exec bash "${SCRIPT_DIR}/install_fuse_common.sh" laptop
