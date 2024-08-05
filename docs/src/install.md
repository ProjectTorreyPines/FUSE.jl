# Setup the FUSE environment

## Julia installation

We highly recommend using the [Juliaup](https://github.com/JuliaLang/juliaup) manager to install Julia
* Mac & Linux: `curl -fsSL https://install.julialang.org | sh`
* Windows: `winget install julia -s msstore`

Once installed, restart your termninal to pick-up the `julia` executable

## FUSE installation

FUSE and related packages are registered at the [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). For installation:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("FUSE")
```

## Test FUSE installation
1. Start your Julia interpreter by typing `julia` at the terminal

1. Try importing the FUSE package

   ```julia
   using FUSE
   ```

1. Run the regression tests (optional, this can take 1h+)
    ```julia
    ] test FUSE
    ```

## Install Jupyter-Lab and add Julia kernel to it
1. You will need to [install `jupyter-lab`](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html) if that's not already available on your system 

1. Install the `IJulia` package by running:

   ```bash
   cd ~/.julia/dev/FUSE
   make IJulia
   ```

   !!! note
       Run `make IJulia` every time a new julia version is installed.
       This will setup the single- and multi-thread julia kernels in Jupyter.

1. Start a new Jupyter-lab session (this should open a web-browser page with Jupyter running)

   ```bash
   jupyter-lab
   ```

# Update Julia version
Juliaup will inform you when a new release of Julia is available. For example:

```
The latest version of Julia in the `release` channel is 1.9.0+0.aarch64.apple.darwin14. You currently have `1.8.5+0.aarch64.apple.darwin14` installed. Run:

juliaup update

to install Julia 1.9.0+0.aarch64.apple.darwin14 and update the `release` channel to that version.
```

To update Julia and make FUSE work under the new environment do as follows:

1. Update Julia
   ```bash
   juliaup update
   ```

1. Start julia and add Revise (this is necessary if Revise is imported in your `~/.julia/config/startup.jl`)
   ```julia
   import Pkg
   Pkg.add("Revise")
   ```

1. Remove all old `Manifest.toml` files in the FUSE and related packages (these files are specific to a given Julia version)
   ```bash
   cd FUSE
   make rm_manifests
   ```

1. Install all FUSE dependencies and Jupyter
   ```bash
   cd FUSE
   make install
   make IJulia
   ```
