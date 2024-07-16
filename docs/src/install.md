# Setup the FUSE environment

## Julia installation
We highly recommend using the [Juliaup](https://github.com/JuliaLang/juliaup) manager to install Julia
* Mac & Linux: `curl -fsSL https://install.julialang.org | sh`
* Windows: `winget install julia -s msstore`

Once installed, restart your termninal to pick-up the `julia` executable

## FUSE installation
1. Clone the FUSE repository under the `~/.julia/dev` folder ():

   ```bash
   mkdir -p ~/.julia/dev
   cd ~/.julia/dev
   git clone git@github.com:ProjectTorreyPines/FUSE.jl.git FUSE
   ```

   !!! note
       Julia repositories end with `.jl` but their install folders do not

   !!! note
       To clone the FUSE repository you will need to [setup your public key on git GitHub](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)

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

# Install GACODE
1. Download and Install Xquartz: https://www.xquartz.org

1. Download and Install Xcode version 14.3.1 or below: https://xcodereleases.com

   This will have git: check with `git --version`.
   Typing this will prompt installation of Xcode command line tools.

1. Download and Install Macports: https://www.macports.org/install.php

   !!! note
       You will need to be on a secure wifi network for the macports installation to succeed

1. Install `Emacs`, `gcc`, `mpich`, `fftw`, `netcdf`:
   ```bash
   sudo port install emacs +x11
   sudo port install gcc12
   sudo port select --set gcc mp-gcc12
   sudo port install mpich-gcc12
   sudo port select --set mpi mpich-gcc12-fortran
   sudo port install openmpi-gcc12
   sudo port install fftw-3
   sudo port install fftw-3-long
   sudo port install fftw-3-single
   sudo port install netcdf
   sudo port install netcdf-fortran
   ```

1. Clone gacode:
   ```bash
   git clone git@github.com:gafusion/gacode.git
   ```

1. Set-up gacode settings in `$HOME/.zshrc`:
   ```bash
   export GACODE_PLATFORM=OSX_MONTEREY
   export GACODE_ROOT=$HOME/gacode
   . ${GACODE_ROOT}/shared/bin/gacode_setup
   ```

   !!! note
       You may need to replace `mpif90-openmpi-mp` with `mpif90-openmpi-gcc12` in `$GACODE_ROOT/platform/build/make.inc.OSX_MONTEREY`, and `mpirun-openmpi-mp` with `mpirun-openmpi-gcc12` in `$GACODE_ROOT/platform/exec/exec.OSX_MONTEREY`

1. For Mac with Apple Silicon:
   ```bash
   conda install -c conda-forge micromamba
   micromamba  install -c smithsp -c conda-forge gacode
   ```

1. Compile:
   ```bash
   cd $GACODE_ROOT ; make
   ```

1. To test that the build is successful, you can run regression tests:
   ```bash
   neo -r
   tglf -r
   cgyro -g reg01
   cgyro -e ./reg01
   ```
