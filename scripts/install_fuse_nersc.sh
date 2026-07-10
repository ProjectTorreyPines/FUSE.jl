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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/install_fuse_common.sh" nersc
