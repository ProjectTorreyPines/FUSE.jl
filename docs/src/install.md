# Setup the FUSE environment

## Julia installation

### Desktop and laptop (juliaup)

We recommend [Juliaup](https://github.com/JuliaLang/juliaup) on personal machines:

* Mac & Linux: `curl -fsSL https://install.julialang.org | sh`
* Windows: `winget install julia -s msstore --accept-source-agreements --accept-package-agreements`

After installation, restart your terminal so the `julia` command is available.

### HPC systems (environment modules)

On many clusters—including **NERSC Perlmutter**—Julia is provided by the site module system instead of juliaup:

```bash
module load julia
julia --version
```

See [On NERSC (Perlmutter)](@ref nersc-install) for depot layout, `fusebot`, and Jupyter notes specific to NERSC.

## FUSE installation {#fuse-installation}

FUSE and related packages are registered at the [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/).
Start Julia (`julia` at the terminal), then:

1. Add the `FuseRegistry` and the `FUSE` package (a fresh install can take 5+ minutes):

   ```julia
   using Pkg
   Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
   Pkg.Registry.add("General")
   Pkg.add("FUSE")
   ```

1. Import FUSE:

   ```julia
   using FUSE
   ```

1. Install the `fusebot` helper (optional but recommended):

   ```julia
   FUSE.install_fusebot()
   ```

   On systems **without** juliaup (typical on HPC), install into a directory already on your `PATH`:

   ```julia
   FUSE.install_fusebot(joinpath(homedir(), ".local", "bin"))
   ```

1. Run the regression tests (optional; can take 1+ hour):

    ```julia
    ] test FUSE
    ```

1. Exit Julia and clone [`FuseExamples`](https://github.com/ProjectTorreyPines/FuseExamples) in your working directory:

   ```bash
   git clone https://github.com/ProjectTorreyPines/FuseExamples
   ```

   Update later with `git fetch && git reset --hard origin/master` (**this discards local changes to those examples**).

## Install Jupyter-Lab with Julia support

### Python on `PATH`

`fusebot install_IJulia` runs `Pkg.build("IJulia")`, which **requires a Python interpreter on your `PATH`**. If build fails with a missing-Python error:

* Activate a conda environment, or
* Install Python and Jupyter, or
* On HPC, `module load python` (site-specific) before running the install.

A known-good optional stack is provided in [`docs/jupyter_environment.yml`](https://github.com/ProjectTorreyPines/FUSE.jl/blob/master/docs/jupyter_environment.yml):

```bash
conda env create -f docs/jupyter_environment.yml
conda activate fuse-jupyter
```

### Install Jupyter / JupyterLab

Install [JupyterLab](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html) if it is not already available.

!!! note "WebIO and Interact"
    The WebIO JupyterLab extension is needed for [`Interact.jl`](https://github.com/JuliaGizmos/Interact.jl?tab=readme-ov-file#usage).

    * JupyterLab 3.x: check with `python -m jupyter labextension list`. You should see `webio-jupyterlab-provider` enabled.
    * Classic Notebook below version 7: check with `python -m jupyter nbextension list` for `webio-jupyter-nbextension`.
    * Notebook **7+** no longer uses classic `nbextension` commands; use the Lab extension only.

    If extensions conflict, pin JupyterLab 3.x (for example `conda install jupyterlab=3.6.7`) and keep WebIO/Interact up to date. Python 3.11 is a good compatibility target for Interact.

### Install IJulia kernels

```bash
fusebot install_IJulia
```

This installs single- and multi-thread Julia kernels. Thread count for the multi-thread kernel follows `JULIA_NUM_THREADS` (default: number of CPUs). Re-run after installing a new Julia version.

Kernels are written directly under `~/.local/share/jupyter/kernels` (or `$JUPYTER_DATA_DIR/kernels`), so registration does **not** depend on the `jupyter` command being on `PATH`. Listing kernels still requires Jupyter:

```bash
python -m jupyter kernelspec list
```

### Start JupyterLab

**Linux and macOS:**

```bash
python -m jupyter lab
```

**Windows** (prefer the Python module form so the correct environment is used):

```powershell
python -m jupyter lab
```

If `python` is not on `PATH`, use the launcher from your conda or Python install, for example `py -m jupyter lab`.

Open the cloned `FuseExamples` folder and run the tutorial notebooks.

## Updating FUSE

1. Watch the [FUSE repository](https://github.com/ProjectTorreyPines/FUSE.jl) for releases.

1. Update like any Julia package:

    ```julia
    ] up
    ```

!!! tip
    See [Managing packages](https://pkgdocs.julialang.org/v1/managing-packages/) in the Julia manual.

## Updating Julia

### With juliaup

1. `juliaup update`
1. In the new Julia version: `using Pkg; Pkg.add("FUSE")`
1. `fusebot install_IJulia`

### With environment modules

Load the new Julia module, reinstall FUSE in that version's depot if needed, then run `fusebot install_IJulia`.

## Cluster-specific notes

```@contents
Pages = ["install_nersc.md", "install_omega.md", "install_saga.md"]
Depth = 2
```
