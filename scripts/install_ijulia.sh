#!/usr/bin/env bash
# Install IJulia kernels and optional WebIO Jupyter extensions.
# Requires Python on PATH so IJulia can locate a Jupyter runtime during `Pkg.build`.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

find_python() {
    for cmd in python3 python; do
        if command -v "$cmd" >/dev/null 2>&1; then
            echo "$cmd"
            return 0
        fi
    done
    return 1
}

if ! PYTHON="$(find_python)"; then
    cat >&2 <<'EOF'
ERROR: No Python interpreter found on PATH.

IJulia needs Python during Pkg.build("IJulia"). Either:
  - activate a conda environment (see docs/jupyter_environment.yml), or
  - load a site Python module (e.g. on NERSC: module load python), then retry.

EOF
    exit 1
fi

export PATH="$(dirname "$(command -v "$PYTHON")"):${PATH}"

if ! "$PYTHON" -c "import sys; print(sys.executable)" >/dev/null 2>&1; then
    echo "ERROR: Python at $(command -v "$PYTHON") is not usable." >&2
    exit 1
fi

echo "Using Python: $($PYTHON -c 'import sys; print(sys.executable)')"

if ! "$PYTHON" -m jupyter --version >/dev/null 2>&1; then
    cat >&2 <<EOF
WARNING: Jupyter is not installed for this Python.
Install a Jupyter stack, for example:

  conda env create -f ${repo_root}/docs/jupyter_environment.yml
  conda activate fuse-jupyter

Kernels will still be written under ~/.local/share/jupyter/kernels (or
JUPYTER_DATA_DIR/kernels) even if the jupyter command is missing.
EOF
fi

julia "${repo_root}/scripts/install_ijulia_kernels.jl"

"$PYTHON" -m pip install --upgrade webio_jupyter_extension || \
    echo "WARNING: could not install webio_jupyter_extension (Interact widgets may not work in Jupyter)."

if "$PYTHON" -m jupyter lab --version >/dev/null 2>&1; then
    "$PYTHON" -m jupyter labextension list || true
fi

NOTEBOOK_MAJOR="$("$PYTHON" -c "
try:
    import notebook
    print(int(notebook.__version__.split('.')[0]))
except Exception:
    print(0)
" 2>/dev/null || echo 0)"

if [[ "${NOTEBOOK_MAJOR}" -ge 1 && "${NOTEBOOK_MAJOR}" -lt 7 ]]; then
    "$PYTHON" -m jupyter nbextension list || true
elif [[ "${NOTEBOOK_MAJOR}" -ge 7 ]]; then
    echo "Notebook ${NOTEBOOK_MAJOR}.x: classic nbextension commands are not used (JupyterLab extensions only)."
else
    echo "notebook package not found; skipping nbextension check (JupyterLab-only installs are fine)."
fi

if "$PYTHON" -m jupyter kernelspec list >/dev/null 2>&1; then
    "$PYTHON" -m jupyter kernelspec list
fi
