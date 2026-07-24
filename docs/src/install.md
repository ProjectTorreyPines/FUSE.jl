# Setup the FUSE environment

This guide walks you through setting up everything you need to run FUSE: the **Julia** language, the **FUSE** package, the **`fusebot`** helper, and an optional **JupyterLab** environment with Julia kernels.

!!! tip "No-install alternative: the FUSE container"
    If you just want to *run* released FUSE (not develop it), a self-contained
    container image with FUSE, all ProjectTorreyPines packages, and a
    precompiled sysimage is published to
    `ghcr.io/projecttorreypines/fuse:<version>`. It starts in seconds on any
    x86_64 machine with Docker/podman
    (`docker run -it ghcr.io/projecttorreypines/fuse:<version>`; on Apple
    Silicon add `--platform linux/amd64`, which runs emulated), and on HPC
    systems via podman-hpc (NERSC) or Singularity (omega). See
    [`deploy/omega-container/README.md`](https://github.com/ProjectTorreyPines/FUSE.jl/blob/master/deploy/omega-container/README.md).

## One-command install

These scripts install FUSE, Revise, fusebot, the Jupyter stack (`fuse` conda env), IJulia kernels, and clone [`FuseExamples`](https://github.com/ProjectTorreyPines/FuseExamples).

They then activate the `fuse` env, run `fusebot install_IJulia` (or `make install_IJulia` / `scripts/install_ijulia.sh` if fusebot fails), and finish by executing the **first three cells** of `FuseExamples/fluxmatcher.ipynb`. A fresh install typically takes **20–40 minutes** (Julia packages + conda + IJulia + the first flux-matcher solve; the notebook cells are often ~6 minutes on one thread).

### Laptop (Linux or macOS)

From any directory on a personal machine. Installs [juliaup](https://github.com/JuliaLang/juliaup) when `julia` is missing and [Miniconda](https://docs.anaconda.com/miniconda/) when `conda` is missing.

```bash
curl -fsSL https://install.julialang.org | sh -s -- -y && \
bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/install_fuse_laptop.sh)
```

From a local `FUSE.jl` clone: `bash scripts/install_fuse_laptop.sh`

Skip the notebook solve with `FUSE_SKIP_VERIFY=1 bash scripts/install_fuse_laptop.sh` if you only want packages + kernels.

### NERSC (Perlmutter)

From a login node (run the depot symlink under [Home quota / memory pressure](@ref nersc-home-quota) first if `$HOME` is under memory pressure). Loads `julia/1.11.7` and `conda` by default, then continues through IJulia, `FuseExamples`, and the fluxmatcher cells (fine on a login node — typically ~6 minutes on one thread for the notebook solve).

```bash
bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/install_fuse_nersc.sh)
```

From a local `FUSE.jl` clone: `bash scripts/install_fuse_nersc.sh`

Override the Julia module with `FUSE_JULIA_MODULE=julia/1.12.0 bash scripts/install_fuse_nersc.sh` when needed. Skip the notebook solve with `FUSE_SKIP_VERIFY=1` if desired. See [On NERSC (Perlmutter)](@ref nersc-install) for depot layout, `fusebot`, and Jupyter notes.

### Windows

From any directory in **PowerShell**. Installs Julia via [winget](https://learn.microsoft.com/en-us/windows/package-manager/winget/) (Microsoft Store) or the Julia App Installer when `julia` is missing, and Miniconda when `conda` is missing. Like the laptop and NERSC scripts, it finishes by running the first three `fluxmatcher.ipynb` cells unless you set `FUSE_SKIP_VERIFY=1`.

```powershell
winget install julia -s msstore --accept-source-agreements --accept-package-agreements --disable-interactivity; `
irm https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/install_fuse_windows.ps1 | iex
```

From a local `FUSE.jl` clone: `.\scripts\install_fuse_windows.ps1`

If `winget` is unavailable, install Julia with `Add-AppxPackage -AppInstallerFile https://install.julialang.org/Julia.appinstaller`, open a new terminal, then run the install script.

### Re-verify `fluxmatcher.ipynb`

The one-command installs already run the **first three cells** of [`FuseExamples/fluxmatcher.ipynb`](https://github.com/ProjectTorreyPines/FuseExamples/blob/master/fluxmatcher.ipynb). To re-run them later:

* **Cell 0** (code): `using Revise`, `using Plots`, `using FUSE`
* **Cell 1** (markdown): flux-matcher introduction (checked for presence, not executed)
* **Cell 2** (code): flux-matches the DIII-D L-mode case — often ~6 minutes on one thread the first time (compilation plus the solve)

**Linux, macOS, and NERSC:**

```bash
bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/verify_fluxmatcher_notebook.sh)
```

Or from a `FUSE.jl` clone: `bash scripts/verify_fluxmatcher_notebook.sh`

**Windows:**

```powershell
irm https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/verify_fluxmatcher_notebook.ps1 | iex
```

Or from a `FUSE.jl` clone: `.\scripts\verify_fluxmatcher_notebook.ps1`

Then start JupyterLab in the directory that contains `FuseExamples/`:

```bash
conda activate fuse
python -m jupyter lab
```

```powershell
conda activate fuse
python -m jupyter lab
```

Open `FuseExamples/fluxmatcher.ipynb` and run cells 0–2 in the notebook UI.

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

1. Add the `FuseRegistry` and the `FUSE` package (a fresh install typically takes **15–30 minutes** to
   download and precompile dependencies):

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
       `Pkg.add("FUSE")` already precompiles most dependencies. The first `using FUSE` and your first
       `FUSE.init(...)` smoke test can still take **5–15 minutes** extra while remaining packages and
       actor code compile. This is normal and happens only once per Julia/FUSE version, not on every
       startup.

1. Install the `fusebot` helper (optional but recommended). `fusebot` is a small command-line tool
   bundled with FUSE; its main job is to install the Julia Jupyter kernels (`fusebot install_IJulia`),
   plus a few related utilities. Run `fusebot --help` for the list of user commands. Install directory
   is picked by `install_fusebot()` automatically - the juliaup `bin` directory on a laptop, or
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
   ini, act = FUSE.case_parameters(:FPP)
   dd = FUSE.init(ini, act)   # if this completes without error, your install is working
   ```

1. Run the regression tests (optional; can take 1+ hour). At the Julia prompt, typing `]` switches to the package prompt:

    ```julia
    ] test FUSE
    ```

1. Exit Julia and clone [`FuseExamples`](https://github.com/ProjectTorreyPines/FuseExamples) in your working directory:

   ```bash
   git clone https://github.com/ProjectTorreyPines/FuseExamples
   ```

   Update later with `git fetch && git reset --hard origin/master` (**this discards local changes to those examples**).

!!! tip "Laptop quick-start (all-in-one)"
    On a personal machine with juliaup, the full sequence is (after typing `julia` at the terminal):

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

A known-good optional stack is provided in [`docs/jupyter_environment.yml`](https://github.com/ProjectTorreyPines/FUSE.jl/blob/master/docs/jupyter_environment.yml).
That file ships inside the FUSE package, so its location depends on whether FUSE is installed under
`~/.julia/packages` (a normal `Pkg.add`) or `~/.julia/dev` (a `Pkg.develop` checkout). Let Julia resolve
the path for you with `pkgdir(FUSE, ...)` so the commands are copy-paste regardless of where FUSE lives.
Run this in your **terminal** (it calls `julia` for you to locate the file - you do not need to start the
Julia prompt yourself). `julia` must be on your `PATH`: with juliaup it already is, while on HPC you need
`module load julia` first. The guarded line below loads the module on HPC and is a harmless no-op on a laptop:

```bash
command -v module >/dev/null && module load julia   # HPC only; skipped automatically on a laptop
conda env create -f "$(julia -e 'using FUSE; print(pkgdir(FUSE, "docs", "jupyter_environment.yml"))')"
conda activate fuse
```

!!! tip "Developing FUSE (or any other package)"
    To edit FUSE itself, run:

    ```julia
    using Pkg
    Pkg.develop("FUSE")   # any package registered in the FuseRegistry or General registry
    ```

    This clones the source into the standard editable location `~/.julia/dev/FUSE` and points your
    environment at it (instead of the read-only versioned copy under `~/.julia/packages`), so local edits
    take effect immediately.

    The same works for **any** package you want to develop - for example a FUSE dependency like `IMAS`
    or `TJLF`, or a third-party package. The package does **not** have to be registered; `Pkg.develop`
    accepts three forms:

    * **By name** (`Pkg.develop("Foo")`) - requires the package to be in a registry you have added
      (General or the [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/)), since
      Julia reads the repo URL from there to clone it into `~/.julia/dev/Foo`.
    * **By URL** (`Pkg.develop(url="https://github.com/Org/Foo.jl")`) - for unregistered packages; Julia
      clones directly from the given URL.
    * **By path** (`Pkg.develop(path="/path/to/Foo")`) - no registration needed; points the environment
      at an existing local checkout (clone it yourself first).

    Run `Pkg.free("Foo")` to stop developing and return to the registered, versioned copy.

### Install Jupyter / JupyterLab

Install [JupyterLab](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html) if it is not already available.

!!! note "WebIO and Interact"
    The WebIO JupyterLab extension is needed for [`Interact.jl`](https://github.com/JuliaGizmos/Interact.jl?tab=readme-ov-file#usage).

    * JupyterLab 3.x: check with `python -m jupyter labextension list`. You should see `webio-jupyterlab-provider` enabled.
    * Classic Notebook below version 7: check with `python -m jupyter nbextension list` for `webio-jupyter-nbextension`.
    * Notebook **7+** no longer uses classic `nbextension` commands; use the Lab extension only.

    If extensions conflict, pin JupyterLab 3.x (for example `conda install jupyterlab=3.6.7`) and keep WebIO/Interact up to date. Python 3.11 is a good compatibility target for Interact.

### Install IJulia kernels

In your **terminal** (with Python on your `PATH`, and `fusebot` installed earlier):

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

1. Update like any Julia package, from the `]` package prompt (type `]` at the Julia prompt):

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

!!! note "Where to run each command"
    Commands in this guide run in one of two places, indicated by the code-block label and the lead-in text:

    * **Terminal** (your shell) - blocks marked `bash`/`powershell`, e.g. `module load julia`, `git ...`,
      `conda ...`, `fusebot ...`, `python -m jupyter ...`.
    * **Julia prompt** - blocks marked `julia`, run *after* you start Julia's interactive session by
      typing `julia` in the terminal. These use `using`, `Pkg`, `FUSE.`, or the `]` package prompt.

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
