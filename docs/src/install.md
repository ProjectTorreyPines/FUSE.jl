# Setup the FUSE environment

## Julia installation

We highly recommend using the [Juliaup](https://github.com/JuliaLang/juliaup) manager to install Julia
* Mac & Linux: `curl -fsSL https://install.julialang.org | sh`
* Windows: `winget install julia -s msstore`

Once installed, restart your termninal to pick-up the `julia` executable.

## FUSE installation

FUSE and related packages are registered at the [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/).
For installation start your Julia interpreter by typing `julia` at the terminal, then:

1. Add the `FuseRegistry` and the `FUSE` package as you would for any other julia package (for a fresh install this can take 20+ mins):

   ```julia
   using Pkg
   Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
   Pkg.Registry.add("General")
   Pkg.add("FUSE")
   ```

1. Now you should be able to import the FUSE package:

   ```julia
   using FUSE
   ```

1. Install the `fusebot` utility to simplify install/updates later on. Now `fusebot` should be a command that you can type anywhere from the terminal.

   ```julia
   FUSE.install_fusebot()
   ```

1. Run the regression tests (optional, this can take 1h+)

    ```julia
    ] test FUSE
    ```

1. Exit julia and clone [`FUSE examples`](https://github.com/ProjectTorreyPines/FuseExamples) in the current working directory. To see/run those `.ipynb` files, you'll need to use Jupyter-Lab or VScode.

   ```bash
   fusebot install_examples
   ```

## Install Jupyter-Lab with Julia support

1. You will need to [install `jupyter-lab`](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html) if that's not already available on your system 

1. Install the `IJulia` package by running:

   ```bash
   fusebot install_IJulia
   ```

   !!! note
       This will setup the single- and multi-thread julia kernels in Jupyter.

       The number of threads of the multi-threaded julia kernels can be set via the `JULIA_NUM_THREADS` environmental variable.

       This needs to be done every time a new version of Julia is installed.

1. Start a new Jupyter-lab session (this should open a web-browser page with Jupyter running)

   ```bash
   jupyter-lab
   ```

1.  Now you can browse the examples in the `FuseExamples` folder that you have cloned, and take a tour of the example Jupyter notebooks there.

## Updating FUSE

1. Get notified of new FUSE releases by "watching" the [FUSE repo on GitHub](https://github.com/ProjectTorreyPines/FUSE.jl)

1. FUSE is [updated like any other Julia package](https://pkgdocs.julialang.org/v1/managing-packages/#updating):

    ```julia
    ] up
    ```

!!! tip
    Become familiar with how [managing Julia packages](https://pkgdocs.julialang.org/v1/managing-packages/) works.

## Updating Julia

1. Use `juliaup update` to install the latest version of Julia

1. Install FUSE to the new version of Julia

   ```julia
   using Pkg
   Pkg.add("FUSE")
   ```

1. Run `fusebot install_IJulia` to install the Kernel for the latest version of Julia in Jupyter-Lab
