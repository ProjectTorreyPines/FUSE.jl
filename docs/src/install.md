# Setup the FUSE environment

## Julia installation

We highly recommend using the [Juliaup](https://github.com/JuliaLang/juliaup) manager to install Julia
* Mac & Linux: `curl -fsSL https://install.julialang.org | sh`
* Windows: `winget install julia -s msstore --accept-source-agreements --accept-package-agreements`

Once installed, restart your terminal to pick-up the `julia` executable.

## FUSE installation

FUSE and related packages are registered at the [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/).
For installation start your Julia interpreter by typing `julia` at the terminal, then:

1. Add the `FuseRegistry` and the `FUSE` package as you would for any other julia package (for a fresh install this can take 5+ mins):

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
   git clone https://github.com/ProjectTorreyPines/FuseExamples
   ```

   This is a git repository that you are in control of. Do a `git fetch && git reset --hard origin/master` to gather the latest updates (**NOTE: this will wipe out any changes you have made to those examples!**)

## Install Jupyter-Lab with Julia support

1. You will need to [install `jupyter-lab`](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html) if that's not already available on your system

   !!! note
       The WebIO jupyter-lab extension is needed for the [`Interact.jl`](https://github.com/JuliaGizmos/Interact.jl?tab=readme-ov-file#usage) package to work.

       Make sure WebIO is working with `jupyter labextension list`. If it is working properly, you should see something like: webio-jupyterlab-provider v0.1.0 enabled OK (python, webio_jupyter_extension). This can also be checked with `jupyter nbextension list`, which should show something like: webio-jupyter-nbextension/nbextension  enabled.

       If the extension has compatibility issues, install an older version of Jupyter (eg. `conda install jupyterlab=3.6.7`). Also ensure that the WebIO and Interact packages are fully up-to-date, and restart the notebook session before testing that it works. Finally, it may be necessary to downgrade your system's version of Python - the recommended version for compatibility with Interact is 3.11.11.

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

## Working on GA clusters

```@contents
Pages = ["install_omega.md", "install_saga.md"]
Depth = 2
```
