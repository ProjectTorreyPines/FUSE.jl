# FUSE.jl

Getting started
===============

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

Key FUSE concepts
=================

In FUSE:
* Data is organized according to the [ITER IMAS hierarchy](https://gafusion.github.io/omas/schema.html)
  * `dd = IMAS.dd()` (which stands for "data dictionary" is the root of the data structure
* To first populate the data dictionary and run any actors FUSE one can:
  * Manually fill the `dd` data structure
  * Read in an existing OMAS JSON data structure
  * Use the `FUSE.Init(dd, ini)` method to populates `dd` starting from 0D `ini` parameters (same spirit of OMFIT's PRO_create module)
* Physics and engineering "actors" are the building blocks of simulations and operate **exclusively** on the `dd data dictionary and their functionality is controlled via `act` parameters
* Both the `ini` and `act` parameters structures can be thought of a glorified namelists, which can be:
  * Filled starting from scratch
  * Taken (and modified) from pre-defined cases (see under the `FUSE/cases/` folder)
  * Taken (and modified) from GASC outputs

![image](https://user-images.githubusercontent.com/1537880/167070559-aeb20212-de01-4fff-ba68-4ebe70cc2b18.png)
