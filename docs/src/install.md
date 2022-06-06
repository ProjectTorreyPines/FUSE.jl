# Setup the FUSE environment

## Julia installation
* Download the latest stable Julia version: [https://julialang.org](https://julialang.org)
* Setup your `$PATH` environment variable
  - Under OSX edit your `~/.zshrc` and add ` PATH=/Applications/Julia-1.7.app/Contents/Resources/julia/bin:$PATH`
  - Under Linux edit your `~/.bashrc` and add ` PATH=/PATH_TO_PARENT_JULIA_FOLDER/julia/bin:$PATH`
* Test that your julia interpreter starts by typing `julia` at the terminal

## Install FUSE packages
1. Clone the FUSE repository under the `~/.julia/dev` folder (note that the repository ends with `.jl` but the install folder does not):
   ```bash
   mkdir -p ~/.julia/dev
   cd ~/.julia/dev
   git clone git@github.com:ProjectTorreyPines/FUSE.jl.git FUSE
   ```

!!! note
    In order to git clone the repository you will need to setup you need to upload your public key to git see [git_documentation](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)

1. Add FUSE and its dependencies to the julia environment (this may take a few minutes):
    ```bash
    cd ~/.julia/dev/FUSE
    make develop
    ```

!!! note
    To pull latest changes for FUSE and all of its dependencies:

    ```bash
    cd ~/.julia/dev/FUSE
    make update
    ```

## Test FUSE installation
1. Start your julia interpreter by typing `julia` at the terminal
1. Try importing the FUSE package
   ```julia
   using FUSE
   ```

## Install Jupyter and add Julia kernel to it
!!! note
    On a fresh install you will need to install python and jupyter notebooks (via miniconda)
1. Install the IJulia package manually or by running:
   ```bash
   cd ~/.julia/dev/FUSE
   make IJulia
   ```
1. Start a Jupyter-lab session (this should open a web-browser page with Jupyter running)
   ```bash
   jupyter-lab
   ```
