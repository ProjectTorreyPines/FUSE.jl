# On NERSC (Perlmutter) {#nersc-install}

These instructions are for **personal** FUSE installs on NERSC login or compute nodes.
For the pre-built `module load fuse` environment maintained by the FUSE team, see the
[Perlmutter deployment scripts](https://github.com/ProjectTorreyPines/FUSE.jl/tree/master/deploy/perlmutter).

## Julia

NERSC provides Julia through the module system. **Do not use juliaup** for a typical NERSC workflow.

```bash
module load julia    # or: module avail julia && module load julia/<version>
julia --version
```

Use a project depot on `$SCRATCH` if your home quota is tight:

```bash
export JULIA_DEPOT_PATH="${SCRATCH}/.julia:${JULIA_DEPOT_PATH:-}"
mkdir -p "${SCRATCH}/.julia"
```

## FUSE

Start Julia and follow the [FUSE installation](@ref fuse-installation) steps (registry + `Pkg.add("FUSE")`).

Install `fusebot` and add it to your shell startup file (the juliaup-style step that `module load julia` does **not** do for you):

```julia
using FUSE
FUSE.install_fusebot(; setup_shell=true)
```

| | **juliaup (laptop)** | **`module load julia` (NERSC)** |
|---|---|---|
| Puts `julia` on `PATH` | Yes, via installer hook in `~/.bashrc` | Yes, **only while the module is loaded** |
| Puts `fusebot` on `PATH` | Yes, if installed next to juliaup | **No** — you must install to `~/.local/bin` (or similar) |
| Persistent user-tool `PATH` | `~/.juliaup/bin` added automatically | Use `FUSE.install_fusebot(; setup_shell=true)` once |

`setup_shell=true` appends a marked block to `~/.bashrc` / `~/.zshrc`:

```bash
# FUSE fusebot PATH
export PATH="$HOME/.local/bin:$PATH"
```

In the **current** shell (before re-login), either run `export PATH="$HOME/.local/bin:$PATH"` or `module load julia` and then use the full path once: `~/.local/bin/fusebot`.

Alternatively, run `FUSE.setup_fusebot_shell!()` if you already installed `fusebot` but it is not on `PATH`.

## Python and Jupyter

`fusebot install_IJulia` (and `Pkg.build("IJulia")`) need a **Python interpreter on `PATH`**.
On NERSC this often means loading a Python module or activating a personal conda environment **before** running the install:

```bash
module load python   # if available on your system
# or:
conda activate fuse-jupyter   # see docs/jupyter_environment.yml
```

Create the optional conda environment from the FUSE repository (clone or use your checkout path):

```bash
conda env create -f /path/to/FUSE.jl/docs/jupyter_environment.yml
conda activate fuse-jupyter
fusebot install_IJulia
```

Kernels are registered under `~/.local/share/jupyter/kernels` (or `$JUPYTER_DATA_DIR/kernels` if set),
without requiring the `jupyter` executable on `PATH`.

## Running Jupyter on Perlmutter

From a login or interactive node, after kernels are installed:

```bash
python -m jupyter lab --no-browser --port 8888
```

Use SSH port forwarding from your laptop to reach the server, as you would for other HPC systems.

## Pre-built FUSE module (optional)

Release builds installed under `/global/common/software/m3739/perlmutter/fuse` are loaded with:

```bash
module use /global/common/software/m3739/perlmutter/fuse/modules/fuse
module load fuse
```

That environment includes sysimage-accelerated `fuse` and `julia` wrappers and pre-registered Jupyter kernels.
