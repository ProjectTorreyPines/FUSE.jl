# Installation

1. Install Julia
   * https://julialang.org/
   * under OSX edit your `~/.zshrc` and add ` PATH=/Applications/Julia-1.7.app/Contents/Resources/julia/bin:$PATH`
   * start your julia interpreter with `julia`

1. Add GA IR&D Julia registry
   ```julia
   ]  # type "]" to enter package mode
   pkg> add Revise
   pkg> registry add git@github.com:ProjectTorreyPines/GAregistry.git
   ```

1. Add FUSE as development package (this will clone FUSE.jl in the the `FUSE` development folders under `~/.julia/dev`)
   ```julia
   pkg> develop FUSE
   ```
    
1. Add FUSE dependencies as development packages:
   ```julia
   ; # type ";" to enter shell mode
   shell> cd ~/.julia/dev/FUSE
   shell> make develop
   ```

   NOTE: In the future, to pull latest changes for FUSE its dependencies:
   ```
   ; # type ";" to enter shell mode
   shell> cd ~/.julia/dev/FUSE
   shell> make update
   ```
    
1. Add Julia kernel to Jupyter-lab:
   ```
   ; # type ";" to enter shell mode
   shell> cd ~/.julia/dev/FUSE
   shell> make IJulia
   ```