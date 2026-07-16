#!/usr/bin/env bash
# One-shot FUSE install for NERSC Perlmutter login nodes.
#
# Copy-paste (from any working directory, e.g. $HOME or $PSCRATCH):
#   bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/install_fuse_nersc.sh)
#
# Or, from a FUSE.jl clone:
#   bash scripts/install_fuse_nersc.sh
#
# Uses module load julia/1.11.7 and module load conda by default.
# If ~/.local/bin is not writable, fusebot is installed under ~/.local/shared/bin
# and underlying make targets are used when fusebot is not on PATH.

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
exec bash "${SCRIPT_DIR}/install_fuse_common.sh" nersc
