# [On NERSC (Perlmutter)](@id nersc-install)

These instructions are for **personal** FUSE installs on NERSC login or compute nodes.
For the pre-built `module load fuse` environment maintained by the FUSE team, see the
[Perlmutter deployment scripts](https://github.com/ProjectTorreyPines/FUSE.jl/tree/master/deploy/perlmutter).

## Julia

NERSC provides Julia through the module system. **Do not use juliaup** for a typical NERSC workflow.

```bash
module load julia    # or: module avail julia && module load julia/<version>
julia --version
```

### Home quota / memory pressure

If your home directory is tight on space or inode quota, keep the full Julia depot off `$HOME` by
symlinking `~/.julia` to `$PSCRATCH`:

!!! warning "PSCRATCH is purged"
    NERSC periodically purges unused files on `$PSCRATCH`. Treat it as scratch storage, not
    long-term archive. Back up anything you need to keep (for example under `$HOME`, `$CFS`, or
    project storage) before it is removed.

!!! tip "Keep the depot from being purged"
    NERSC resets the purge timer whenever a file is accessed, so periodically refreshing the access
    time on the depot keeps it alive (sufficient if you log in at least about every two months):

    ```bash
    find $PSCRATCH/.julia -exec touch -a {} +   # reset atime so the depot is not purged
    ```

    This is a mitigation, not a guarantee - still back up anything you cannot afford to lose.

```bash
# If ~/.julia already exists as a directory, move it first
if [ -d ~/.julia ] && [ ! -L ~/.julia ]; then
    mv ~/.julia $PSCRATCH/.julia
else
    mkdir -p $PSCRATCH/.julia
fi

# Create the symlink
ln -s $PSCRATCH/.julia ~/.julia

# Verify
ls -la ~/.julia
```

After this, Julia and `Pkg` use `$PSCRATCH/.julia` transparently. You do not need a separate
`JULIA_DEPOT_PATH` unless you want multiple depots (for example a read-only site depot).

Alternatively, without moving `~/.julia`, you can prepend a scratch depot:

```bash
export JULIA_DEPOT_PATH="${SCRATCH}/.julia:${JULIA_DEPOT_PATH:-}"
mkdir -p "${SCRATCH}/.julia"
```

## FUSE

Start Julia and follow the [FUSE installation](@ref fuse-installation) steps (registry + `Pkg.add("FUSE")`).

Then install `fusebot`. You do **not** pick the directory: under `module load julia` it is installed
automatically to `~/.local/bin`, and `setup_shell=true` adds that directory to your shell `PATH`:

```julia
using FUSE
FUSE.install_fusebot(; setup_shell=true)   # installs to ~/.local/bin and adds it to your PATH
```

!!! note "Why `setup_shell=true` is needed on NERSC"
    `module load julia` puts `julia` on your `PATH` **only while the module is loaded** and never adds
    user tools like `fusebot` (unlike juliaup, which configures `~/.juliaup/bin` automatically). Running
    this once writes a marked block to `~/.bashrc` / `~/.zshrc` so `fusebot` is on `PATH` in new shells:

    ```bash
    # FUSE fusebot PATH
    export PATH="$HOME/.local/bin:$PATH"
    ```

!!! tip "Why `~/.local/bin`?"
    `~/.local/bin` lives in `$HOME` (global homes), which is mounted on both login and compute nodes and,
    unlike `$PSCRATCH`, is persistent and not purged (it does have a quota). Run the install from a login
    node. To use a different location, pass it explicitly: `FUSE.install_fusebot("/custom/bin"; setup_shell=true)`.

To use `fusebot` in the **current** shell before re-login, run `export PATH="$HOME/.local/bin:$PATH"`
(or call it once by full path: `~/.local/bin/fusebot`). If you already installed `fusebot` but it is not
on `PATH`, run `FUSE.setup_fusebot_shell!()`.

!!! tip "NERSC quick-start (all-in-one)"
    From a login node, after the one-time `~/.julia` depot setup above:

    ```bash
    module load julia
    julia
    ```

    then, inside Julia:

    ```julia
    using Pkg
    Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
    Pkg.Registry.add("General")
    Pkg.add("FUSE")
    using FUSE
    FUSE.install_fusebot(; setup_shell=true)
    ```

## Python and Jupyter

`fusebot install_IJulia` (and `Pkg.build("IJulia")`) need a **Python interpreter on `PATH`**.
On NERSC this often means loading a Python module or activating a personal conda environment **before** running the install:

```bash
module load python   # if available on your system
# or:
conda activate fuse   # see docs/jupyter_environment.yml
```

Create the optional conda environment from the FUSE repository (clone or use your checkout path):

```bash
conda env create -f /path/to/FUSE.jl/docs/jupyter_environment.yml
conda activate fuse
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
