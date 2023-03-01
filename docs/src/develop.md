# Development

## Overview

All FUSE-related development occurs under the [GitHub ProjectTorreyPines organization](https://github.com/ProjectTorreyPines), which lives under the GeneralAtomics GitHub enterprise account.

### Tracking progress

We use a GitHub to [track progress with FUSE developments](https://github.com/orgs/ProjectTorreyPines/projects/2/views/1) and have a birds-eye view across the different repositories used by the FUSE project.

### Packages organization

The FUSE project is built upon different Julia packages. Several of these are managed by GA-MFE, and they all reside in the [https://github.com/ProjectTorreyPines](https://github.com/ProjectTorreyPines) repository.

## How to add/modify entries in `dd`

1. add/edit Json files in the `IMASDD/data_structures_extra` folder
2. run `IMASDD/src/generate_dd.jl`

!!! note
    The `dd` data structure is defined as a Julia `struct`. Like all `struct` re-definitions, changes to the `dd` data structure will requires your Julia interpreters to be restarted to pick-up the updates.

## How to add/modify entries in `ini` and `act`

The functinoality of the `ini` and `act` parameters is implemented in the [SimulationParameters.jl](https://github.com/ProjectTorreyPines/SimulationParameters.jl) package.

* The `ini` parameters are all defined in the `FUSE/src/parameters_init.jl` file. Add/edit entries there.
* The `act` parameters of each actor are defined where that actor is defined. Add/edit entries there.

## How to write a new actor

Actors are grouped in four main abstract types:

```julia
abstract type FacilityAbstractActor <: AbstractActor end
abstract type ReactorAbstractActor <: AbstractActor end
abstract type HCDAbstractActor <: AbstractActor end
abstract type PlasmaAbstractActor <: AbstractActor end
```

The definition of each FUSE actor follows a well defined pattern.
**DO NOT** deviate from this pattern. This is important to ensure modularity and compostability of the actors.

```julia
# Defintion of the actor structure
mutable struct ActorNAME <: ???AbstractActor
    dd::IMAS.dd
    par::ParametersActor  # Actors must carry with them the parameters they are run with
    something_else::??? # Some actors may want to carry something else with them
    # Inner constructor for the actor starting from `dd` and `par` (we generally refer to `par` as `act.ActorNAME`)
    function ActorNAME(dd::IMAS.dd, par::FUSEparameters__ActorNAME; kw...)
        logging_actor_init(ActorNAME)
        par = par(kw...)
        return new(dd, par, something_else)
    end
end

# Definition of the `act` parameters relevant to the actor
# NOTE: To create a `ActorNAME` in `act` you'll have to add these to the FUSE/src/parameters_actors.jl file
Base.@kwdef mutable struct FUSEparameters__ActorNAME{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    length::Entry{T} = Entry(T, "m", "Some decription") # it's ok not to have a default, it forces users to think about what a parameter should be
    verbose::Entry{Bool} = Entry(Bool, "", "Some other decription"; default=true)
    switch::Switch{Symbol} = Switch(Symbol, [:option_a, :option_b], "", "user can only select one of these"; default=:option_a)
end

# Constructor with with `dd` and `act` as arguments will actually run the actor!
# That's how users will mostly run this actor.
# This does not change, and it's always the same for all actors
"""
    ActorNAME(dd::IMAS.dd, act::ParametersAllActors; kw...)

What does this actor do...
"""
function ActorNAME(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorNAME(kw...) # this makes a local copy of `act.ActorNAME` and overrides it with keywords that the user may have passed
    actor = ActorNAME(dd, par) # instantiate the actor (see function below)
    step(actor)                # run the actor
    finalize(actor)            # finalize
    return actor
end

# define `_step` function for this actor (this is where most of the action occurs)
# note the leading underscore (use the `_step()` and not `step()` for the FUSE logging system to work with your actor)
# `_step()` should not take any argument besides the actor itself
function _step(actor::ActorNAME)
    ...
    return actor # _step() should always return the actor
end

# define `_finalize` function for this actor (this is where typically data gets written in `dd` if that does happen already at the `step`)
# note the leading underscore (use the `_finalize()` and not `finalize()` for the FUSE logging system to work with your actor)
# `_finalize()` should not take any argument besides the actor itself
function _finalize(actor::ActorNAME)
    ...
    return actor # _finalize() should always return the actor
end
```

## How to build the documentation

1. To build the documentation, in the `FUSE/docs` folder, start Julia then:
   ```julia
   include("make.jl")
   ```
   !!! note Interactive documentation build
       One can call `include("make.jl")` over and over within the same Julia session to avoid dealing with startup time.

1. Check page by opening `FUSE/docs/build/index.html` page in web-browser.

1. To publish online in the `FUSE` folder:
   ```bash
   make web
   ```

## Tips and more

### Revise.jl
Use [Revise.jl](https://github.com/timholy/Revise.jl) to modify code and use the changes without restarting Julia.
We recommend adding `import Revise` to your `~/.julia/config/startup.jl`.

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
