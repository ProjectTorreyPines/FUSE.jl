# FUSE.jl

Getting started
===============

1. Add GA IR&D Julia registry
    ```julia
    pkg> registry add git@github.com:ProjectTorreyPines/GAregistry.git
    ```

2. Add FUSE as development package (this will clone FUSE.jl in the the `FUSE` development folders under `~/.julia/dev`)
    ```julia
    pkg> develop FUSE
    ```
    
3. Add FUSE dependencies as development packages:
    ```bash
    cd ~/.julia/dev/FUSE
    make develop
    ```

3. Pull latest changes for FUSE its dependencies
    ```bash
    cd ~/.julia/dev/FUSE
    make update
    ```
