#!/usr/bin/env bash
# Step 2 after install: verify fluxmatcher.ipynb cells 0–1 from the command line.
#
# Laptop:
#   bash scripts/verify_fluxmatcher_notebook.sh
#
# NERSC:
#   module load julia/1.11.7 && bash scripts/verify_fluxmatcher_notebook.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

log() { echo "[verify_fluxmatcher] $*"; }

if command -v module >/dev/null 2>&1 && ! command -v julia >/dev/null 2>&1; then
    module load "${FUSE_JULIA_MODULE:-julia/1.11.7}" 2>/dev/null || module load julia
fi

if [[ -f "${HOME}/.juliaup/bin/julia" ]]; then
    export PATH="${HOME}/.juliaup/bin:${PATH}"
fi

command -v julia >/dev/null 2>&1 || { echo "ERROR: julia not on PATH" >&2; exit 1; }

WORK_DIR="${FUSE_WORK_DIR:-${PWD}}"
if [[ ! -f "${WORK_DIR}/FuseExamples/fluxmatcher.ipynb" ]]; then
    echo "ERROR: ${WORK_DIR}/FuseExamples/fluxmatcher.ipynb not found." >&2
    echo "Run the install script first, or set FUSE_WORK_DIR to your FuseExamples parent directory." >&2
    exit 1
fi

cd "${WORK_DIR}"

if [[ -x "${HOME}/.local/miniconda3/envs/fuse/bin/python" ]]; then
    PYTHON="${HOME}/.local/miniconda3/envs/fuse/bin/python"
elif command -v python >/dev/null 2>&1; then
    PYTHON="$(command -v python)"
else
    echo "ERROR: python not on PATH (activate the fuse conda env first)" >&2
    exit 1
fi

"${PYTHON}" - <<'PY'
import json
import sys
from pathlib import Path

path = Path("FuseExamples/fluxmatcher.ipynb")
nb = json.loads(path.read_text())
if len(nb["cells"]) < 2:
    sys.exit("fluxmatcher.ipynb must have at least 2 cells")
if nb["cells"][0]["cell_type"] != "code":
    sys.exit("Cell 0 must be a code cell")
if nb["cells"][1]["cell_type"] != "markdown":
    sys.exit("Cell 1 must be a markdown cell")
text = "".join(nb["cells"][1].get("source", [])).lower()
if "flux-matcher" not in text and "flux matcher" not in text:
    sys.exit("Cell 1 markdown does not mention flux-matcher")
print("[verify_fluxmatcher] Cell 1 (markdown) present")
PY

julia "${SCRIPT_DIR}/verify_fluxmatcher_notebook.jl"
log "fluxmatcher.ipynb cells 0–1 verification PASSED"
