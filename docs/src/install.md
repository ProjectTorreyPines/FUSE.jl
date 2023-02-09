# Setup the FUSE environment

## Julia installation
We highly recommend using the [Juliaup](https://github.com/JuliaLang/juliaup) manager to install Julia
* Mac & Linux: `curl -fsSL https://install.julialang.org | sh`
* Windows: `winget install julia -s msstore`

## Install FUSE packages
1. Clone the FUSE repository under the `~/.julia/dev` folder (note that the repository ends with `.jl` but the install folder does not):

   ```bash
   mkdir -p ~/.julia/dev
   cd ~/.julia/dev
   git clone git@github.com:ProjectTorreyPines/FUSE.jl.git FUSE
   ```

   !!! note
       To clone the FUSE repository you will need to [setup your public key on git GitHub](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)

   !!! note
       If you are installing from an old install and the package update isn't able to resolve itself feel free to remove every folder except the `dev` folder in your `.julia` folder

1. Add FUSE and its dependencies to the julia environment (this may take a few minutes):

   ```bash
   cd ~/.julia/dev/FUSE
   make install
   ```

   !!! note
       To pull latest changes for FUSE and all of its dependencies:

       ```bash
       cd ~/.julia/dev/FUSE
       make update
       ```

## Test FUSE installation
1. Start your Julia interpreter by typing `julia` at the terminal

1. Try importing the FUSE package

   ```julia
   using FUSE
   ```

## Install Jupyter and add Julia kernel to it
1. You will need to [install `jupyter-lab`](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html) if that's not already available on your system 

1. Install the `IJulia` package by running:

   ```bash
   cd ~/.julia/dev/FUSE
   make IJulia
   ```

1. Start a new Jupyter-lab session (this should open a web-browser page with Jupyter running)

   ```bash
   jupyter-lab
   ```

## Troubleshooting
When using the ProjectTorreyPines private Julia [GA registry](https://github.com/ProjectTorreyPines/GAregistry), one may get `SSH host verification` errors when installing and updating Julia packages:

```
SSH host verification: the identity of the server `github.com:22` does not match its known hosts record. Someone could be trying to man-in-the-middle your connection. It is also possible that the server has changed its key, in which case you should check with the server administrator and if they confirm that the key has been changed, update your known hosts file.
```

this can be resolved by telling Julia to use the git command line interface, via the environmental variable:

```julia
export JULIA_PKG_USE_CLI_GIT=true
```

## Install TGLF from Fortran source (optional)


