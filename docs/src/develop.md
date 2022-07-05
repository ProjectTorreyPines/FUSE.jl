# Development

## Overview

All FUSE-related development occurs under the [GitHub ProjectTorreyPines organization](https://github.com/ProjectTorreyPines), which lives under the GeneralAtomics GitHub enterprise account.

### Tracking progress

We use a GitHub to [track progress with FUSE developments](https://github.com/orgs/ProjectTorreyPines/projects/2/views/1) and have a birds-eye view across the different repositories used by the FUSE project.

### Packages organization

The FUSE project is built upon different Julia packages. Several of these are managed by GA-MFE, and they all reside in the [https://github.com/ProjectTorreyPines](https://github.com/ProjectTorreyPines) repository. Of all the packages that are there, three are at the fundation of the FUSE framework itself: `FUSE.jl`, `IMAS.jl`, `CoordinateConventions.jl`, `IMASDD.jl`. All other packages in the ProjectTorreyPines organization fundamentally add physics and technology capabilities in FUSE (typically as actors).

* **FUSE.jl**
  * Actors
  * `ini` parameters
  * `act` parameters
  * Workflows
  * Physics functions
  * Technology functions
  * Utility functions

* **IMAS.jl**
  * Physics functions
  * Physics constants
  * Plotting functions
  * Math functions
  * `dd` expressions

* **CoordinateConventions.jl**
  * Coordinate conventions (COCOS) and transformations

* **IMASDD.jl**
  * `dd` data structure
  * Base `dd` functionality
  * Representation of `dd`
  * Loading/saving of `dd`
  * Utility and math (eg. interpolation)

## Editing stuff

### Add/modify entries in `dd`

1. add/edit Json files in the `IMASDD/data_structures_extra` folder
2. run `IMASDD/src/generate_dd.jl`

!!! note
    The `dd` data structure is defined as a Julia `struct`. Like all `struct` re-definitions, changes to the `dd` data structure will requires your Julia interpreters to be restarted to pick-up the updates.

### Add/modify entries in `ini`

`ini` parameters are of type `ParametersInit`.

1.  The `ParametersInit` are defined all `FUSE/src/parameters_init.jl` file. Add/edit entries there.

### Add/modify entries in `act`

`act` parameters are of type `ParametersActor`.

1. The `ParametersActor` of each actor are defined where that actor is defined. Add/edit entries there.

### What constitutes an Actor

Actors are grouped in four main abstract types:

```julia
abstract type FacilityAbstractActor <: AbstractActor end
abstract type ReactorAbstractActor <: AbstractActor end
abstract type HCDAbstractActor <: AbstractActor end
abstract type PlasmaAbstractActor <: AbstractActor end
```

The definition of each FUSE actor follows a well defined pattern:

```julia
# Defintion of the actor structure
Base.@kwdef mutable struct ActorNAME <: ???AbstractActor
    dd::IMAS.dd
    ...
end

# Definition of the `act` parameters relevant to the actor
function ParametersActor(::Type{Val{:ActorNAME}})
    par = ParametersActor(nothing)
    par.ngrid = Entry(Integer, "", "Grid ..."; default=129)
    ...
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end

# run actor with `dd` and `act` as arguments
"""
    ActorNAME(dd::IMAS.dd, act::ParametersActor; kw...)

What does this actor do...
"""
function ActorNAME(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorNAME(kw...)
    actor = ActorNAME(;dd, par...)
    step(actor)
    finalize(actor)
    return actor
end

# define `step` function for this actor
function step(actor::ActorNAME; ...)
    ...
end

# define `finalize` function for this actor
function finalize(actor::ActorNAME; ...)
    ...
end
```

## Build documentation

1. To build missing example notebooks pages:
   ```bash
   make examples
   ```

1. To build the documentation, in the `FUSE/docs` folder:
   ```julia
   include("make.jl")
   ```
   !!! $note Interactive documentation build
       One can call `include("make.jl")` over and over within the same Julia session to avoid dealing with startup time.

1. Check page by opening `FUSE/docs/build/index.html` page in web-browser.

1. To publish online
   ```bash
   make web
   ```

## Tips and more

### Development in VSCode

[VSCode](https://code.visualstudio.com) is an excellent development environment for Julia.

FUSE uses the following VSCode settings for formatting the Julia code:
```json
{
    "files.autoSave": "onFocusChange",
    "workbench.tree.indent": 24,
    "[julia]": {
        "editor.defaultFormatter": "singularitti.vscode-julia-formatter"
    },
    "juliaFormatter.margin": 160,
    "juliaFormatter.alwaysForIn": true,
    "juliaFormatter.annotateUntypedFieldsWithAny": false,
    "juliaFormatter.whitespaceInKwargs": false,
    "juliaFormatter.overwriteFlags": true,
    "juliaFormatter.alwaysUseReturn": true,
}
```

!!! note
    To add these settings to VSCode add these lines to: `<Command> + <shift> + p` -> `Preferences: Open Settings (JSON)`

In addition, we suggest enabling the VSCode `autoSave` feature and increasing the `workbench tree indent` to 24

To format Julia you will need to install `Julia Language Support` under the extensions tab (`<Command> + <shift> + x`)

!!! tip
    To see what is precompiled at runtime, you can add a Julia kernel with the `trace-compile` option to Jupyter
    ```julia
    import IJulia
    IJulia.installkernel("Julia tracecompile", "--trace-compile=stderr")
    ```
    Then select the `Julia tracecompile` in jupyter-lab

!!! note
    When pushing a jupyter notebook make sure that the output is cleared 
