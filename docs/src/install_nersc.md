# [On NERSC (Perlmutter)](@id nersc-install)

These instructions are for **personal** FUSE installs on NERSC login or compute nodes.
For the pre-built `module load fuse` environment maintained by the FUSE team, see the
[Perlmutter deployment scripts](https://github.com/ProjectTorreyPines/FUSE.jl/tree/master/deploy/perlmutter).

## NERSC one-command install

From a login node (run the depot symlink under [Home quota / memory pressure](@ref nersc-home-quota) first if `$HOME` is tight). Loads `julia/1.11.7` and `conda`, installs FUSE, Revise, fusebot (under `~/.local/bin` or `~/.local/shared/bin` if `bin` is not writable), the Jupyter conda environment, IJulia kernels, and clones [`FuseExamples`](https://github.com/ProjectTorreyPines/FuseExamples). All steps run from the shell — no Julia REPL.

```bash
bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/install_fuse_nersc.sh)
```

From a local `FUSE.jl` clone: `bash scripts/install_fuse_nersc.sh`

Override the Julia module with `FUSE_JULIA_MODULE=julia/1.12.0 bash scripts/install_fuse_nersc.sh` when needed.

If `fusebot` cannot be installed (for example when `~/.local/bin` does not exist), the script runs the same underlying `make` targets (`install_IJulia`, `install_examples`, …) directly from the FUSE package directory.

### Step 2: verify `fluxmatcher.ipynb`

Confirms **cell 0** (`using Revise`, `using Plots`, `using FUSE`) and **cell 1** (markdown) in `FuseExamples/fluxmatcher.ipynb`:

```bash
module load julia/1.11.7
bash <(curl -fsSL https://raw.githubusercontent.com/ProjectTorreyPines/FUSE.jl/master/scripts/verify_fluxmatcher_notebook.sh)
```

Or from a `FUSE.jl` clone: `bash scripts/verify_fluxmatcher_notebook.sh`

Then start JupyterLab (with `module load conda` and `conda activate fuse`) and open `FuseExamples/fluxmatcher.ipynb`.

## Julia

NERSC provides Julia through the module system. **Do not use juliaup** for a typical NERSC workflow.

```bash
module load julia    # or: module avail julia && module load julia/<version>
julia --version
```

If your home directory is under memory pressure, move the Julia depot off `$HOME` onto
`$PSCRATCH` - see [Home quota / memory pressure](@ref nersc-home-quota) under Troubleshooting below.

## FUSE

Start Julia and follow the [FUSE installation](@ref fuse-installation) steps (registry + `Pkg.add("FUSE")`).

Then install `fusebot`. You do **not** pick the directory: under `module load julia` it is installed
automatically to `~/.local/bin`, and `setup_shell=true` adds that directory to your shell `PATH`:

```julia
using FUSE
FUSE.install_fusebot(; setup_shell=true)   # installs to ~/.local/bin and adds it to your PATH
```

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

## Python and Jupyter

`fusebot install_IJulia` (and `Pkg.build("IJulia")`) need a **Python interpreter on `PATH`**. The simplest
option is the optional conda environment bundled with FUSE, which provides Python and JupyterLab. Let Julia
locate the bundled `jupyter_environment.yml` so the command is copy-paste regardless of where FUSE is
installed (`module load julia` makes `julia` available):

```bash
module load julia
conda env create -f "$(julia -e 'using FUSE; print(pkgdir(FUSE, "docs", "jupyter_environment.yml"))')"
conda activate fuse
fusebot install_IJulia
```

Any Python on `PATH` works if you prefer not to use this environment - for example `module load python`
or activating your own conda environment before running `fusebot install_IJulia`.

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

## Troubleshooting

### [Home quota / memory pressure](@id nersc-home-quota)

If your home directory is tight on space or inode quota, keep the full Julia depot off `$HOME` by
symlinking `~/.julia` to `$PSCRATCH`:

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

    To automate this, add the following to your `~/.bashrc`. It runs only on interactive shells, at
    most once per day, and in the background so it never delays your prompt:

    ```bash
    # Refresh atime on the Julia depot so the PSCRATCH purge does not remove it.
    if [[ $- == *i* && -d "$PSCRATCH/.julia" ]]; then
        _depot_stamp="$PSCRATCH/.julia/.last_atime_refresh"
        if [[ ! -e "$_depot_stamp" || -n $(find "$_depot_stamp" -mtime +1 -print 2>/dev/null) ]]; then
            ( touch "$_depot_stamp"
              find "$PSCRATCH/.julia" -exec touch -a {} + ) >/dev/null 2>&1 &
        fi
        unset _depot_stamp
    fi
    ```

    This is a mitigation, not a guarantee - still back up anything you cannot afford to lose, and note
    that if you go longer than the full purge window without logging in, the depot can still be removed.
