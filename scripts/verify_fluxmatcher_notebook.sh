#!/usr/bin/env bash
# Verify fluxmatcher.ipynb cells 0–2 from the command line.
# Also the final step of the laptop and NERSC one-command installs.
#
# Laptop:
#   bash scripts/verify_fluxmatcher_notebook.sh
#
# NERSC:
#   module load julia/1.11.7 && bash scripts/verify_fluxmatcher_notebook.sh

set -euo pipefail

SCRIPT_BASE_URL="${FUSE_SCRIPT_BASE_URL:-https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts}"

resolve_script_dir() {
    local self_dir
    self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -f "${self_dir}/verify_fluxmatcher_notebook.jl" ]]; then
        echo "${self_dir}"
        return 0
    fi

    # curl | bash / process-substitution: download the Julia runner next to us.
    local bundle_dir="${TMPDIR:-/tmp}/fuse-verify-$$"
    mkdir -p "${bundle_dir}"
    curl -fsSL "${SCRIPT_BASE_URL}/verify_fluxmatcher_notebook.jl" \
        -o "${bundle_dir}/verify_fluxmatcher_notebook.jl"
    echo "${bundle_dir}"
}

SCRIPT_DIR="$(resolve_script_dir)"

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

# Extract the code from the first three cells (cells 0 & 2; cell 1 is markdown)
# into a temp file that verify_fluxmatcher_notebook.jl runs, so the test stays
# faithful to whatever the notebook actually contains.
export FUSE_VERIFY_CELLS="$(mktemp)"
trap 'rm -f "${FUSE_VERIFY_CELLS}"' EXIT

"${PYTHON}" - <<'PY'
import json
import os
import sys
from pathlib import Path

out_path = Path(os.environ["FUSE_VERIFY_CELLS"])
path = Path("FuseExamples/fluxmatcher.ipynb")
nb = json.loads(path.read_text())
cells = nb["cells"]
if len(cells) < 3:
    sys.exit("fluxmatcher.ipynb must have at least 3 cells")
if cells[0]["cell_type"] != "code":
    sys.exit("Cell 0 must be a code cell")
if cells[1]["cell_type"] != "markdown":
    sys.exit("Cell 1 must be a markdown cell")
if cells[2]["cell_type"] != "code":
    sys.exit("Cell 2 must be a code cell")
text = "".join(cells[1].get("source", [])).lower()
if "flux-matcher" not in text and "flux matcher" not in text:
    sys.exit("Cell 1 markdown does not mention flux-matcher")

chunks = []
for i, cell in enumerate(cells[:3]):
    if cell["cell_type"] != "code":
        continue
    chunks.append("# --- fluxmatcher.ipynb cell %d ---" % i)
    chunks.append("".join(cell.get("source", [])))
out_path.write_text("\n".join(chunks) + "\n")
print("[verify_fluxmatcher] Cells 0-2 present (cell 1 markdown, cells 0 & 2 code)")
PY

julia "${SCRIPT_DIR}/verify_fluxmatcher_notebook.jl"
log "fluxmatcher.ipynb cells 0–2 verification PASSED"
