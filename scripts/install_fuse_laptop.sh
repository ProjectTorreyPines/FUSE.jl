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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/install_fuse_common.sh" laptop
