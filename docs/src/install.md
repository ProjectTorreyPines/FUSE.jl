# Setup the FUSE environment

This guide walks you through setting up everything you need to run FUSE: the **Julia** language, the **FUSE** package, the **`fusebot`** helper, and an optional **JupyterLab** environment with Julia kernels.

!!! note "Prerequisites"
    * **Julia 1.10 or newer** (installed below via juliaup or your site's module system).
    * **`git`** on your `PATH` (to clone the examples).
    * Works on **Linux, macOS, and Windows**, as well as HPC systems.
    * Budget **~15 minutes** for the install, plus a one-time precompilation the first time you run `using FUSE`.

!!! tip "What is `fusebot`?"
    `fusebot` is a small command-line helper bundled with FUSE. Its main job is to install the Julia
    Jupyter kernels (`fusebot install_IJulia`), plus a few related utilities. It is optional but recommended.

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

## [FUSE installation](@id fuse-installation)

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

   !!! note "First import is slow"
       The first `using FUSE` (and your first simulation run) triggers precompilation that can take
       several minutes. This is normal and happens only once per Julia/FUSE version, not on every startup.

1. Install the `fusebot` helper (optional but recommended). You do **not** choose the directory:
   `install_fusebot()` picks it automatically - the juliaup `bin` directory on a laptop, or
   `~/.local/bin` under `module load julia` on HPC:

   ```julia
   FUSE.install_fusebot()                                  # auto: juliaup bin, or ~/.local/bin on HPC
   FUSE.install_fusebot(; setup_shell=true)                # HPC: also add ~/.local/bin to your shell PATH
   FUSE.install_fusebot("/custom/bin"; setup_shell=true)   # optional: explicit install directory
   ```

   On HPC (`module load julia`), the site Julia module does not put user tools on `PATH` the way juliaup
   does, so pass `setup_shell=true` once to add the install directory to your shell startup file.
   See [On NERSC (Perlmutter)](@ref nersc-install) for details.

1. Verify the install works with a quick smoke test:

   ```julia
   using FUSE
   ini, act = FUSE.case_parameters(:ITER)
   dd = FUSE.init(ini, act)   # if this completes without error, your install is working
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

!!! tip "Laptop quick-start (all-in-one)"
    On a personal machine with juliaup, the full sequence is:

    ```julia
    using Pkg
    Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
    Pkg.Registry.add("General")
    Pkg.add("FUSE")
    using FUSE
    FUSE.install_fusebot()
    ```

## Install Jupyter-Lab with Julia support

### Python on `PATH`

`fusebot install_IJulia` runs `Pkg.build("IJulia")`, which **requires a Python interpreter on your `PATH`**. If build fails with a missing-Python error:

* Activate a conda environment, or
* Install Python and Jupyter, or
* On HPC, `module load python` (site-specific) before running the install.

A known-good optional stack is provided in [`docs/jupyter_environment.yml`](https://github.com/ProjectTorreyPines/FUSE.jl/blob/master/docs/jupyter_environment.yml):

```bash
conda env create -f docs/jupyter_environment.yml
conda activate fuse
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

## Troubleshooting

!!! warning "`fusebot: command not found`"
    The install directory is not on your `PATH` in the current shell. Either open a new login shell, or
    add it now and re-run, for example `export PATH="$HOME/.local/bin:$PATH"`. To make this permanent,
    run `FUSE.install_fusebot(; setup_shell=true)` (or `FUSE.setup_fusebot_shell!()` if `fusebot` is
    already installed).

!!! warning "`fusebot install_IJulia` fails with a missing-Python error"
    `Pkg.build("IJulia")` needs a Python interpreter on your `PATH`. Activate a conda environment,
    install Python/Jupyter, or on HPC `module load python` (site-specific) before re-running.

!!! warning "`Pkg.Registry.add` fails"
    Make sure `git` is installed and that you can reach GitHub. Behind a proxy or offline node, configure
    your proxy first. You can re-run the registry/add commands; they are safe to repeat.

!!! warning "Wrong Jupyter kernel"
    `fusebot install_IJulia` registers single- and multi-thread Julia kernels. In JupyterLab pick the
    kernel matching the Julia version you installed FUSE into; list them with
    `python -m jupyter kernelspec list`.

## Next steps

* Follow the [introductory tutorial](https://fuse.help/dev/tutorial.html).
* Explore the [`FuseExamples`](https://github.com/ProjectTorreyPines/FuseExamples) notebooks.
* Stuck or have questions? Join the [Discord community](https://discord.gg/CbjpZH9SKM).
