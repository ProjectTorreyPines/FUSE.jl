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

## How to write IMAS physics functions

IMAS physics functions are structured in IMAS.jl under IMAS/src/physics
All of these functions use and or modify the datastructure (dd) in some way and are used to calculate certain quantities or fill the data structure.

Let's say we want to create a function that calculates the DT fusion and then fill in the core_sources dd with the alpha heating from that source. 

Here is an example of writing it in a good way:

```julia
function DT_fusion_source!(dd::IMAS.dd)
    return DT_fusion_source!(dd.core_sources, dd.core_profiles)
end

"""
    DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)

Calculates DT fusion heating with an estimation of the alpha slowing down to the ions and electrons, modifies dd.core_sources
"""
function DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)
    cp1d = cp.profiles_1d[]

    polarized_fuel_fraction = getproperty(cp.global_quantities, :polarized_fuel_fraction, 0.0)
    α = alpha_heating(cp1d; polarized_fuel_fraction)
    if sum(α) == 0
        deleteat!(cs.source, "identifier.index" => 6)
        return cs
    end
    ion_to_electron_fraction = sivukhin_fraction(cp1d, 3.5e6, 4.0)

    source = resize!(cs.source, :fusion; allow_multiple_matches=true)
    new_source(
        source,
        source.identifier.index,
        "α",
        cp1d.grid.rho_tor_norm,
        cp1d.grid.volume,
        cp1d.grid.area;
        electrons_energy=α .* (1.0 .- ion_to_electron_fraction),
        total_ion_energy=α .* ion_to_electron_fraction
    )
    return source
end
```

!!! convention
    The documentation string is added to the specialized function `DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)` and the dispatch function `DT_fusion_source!(dd::IMAS.dd)` is added on top of the function


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
# Definition of the `act` parameters relevant to the actor
Base.@kwdef mutable struct FUSEparameters__ActorNAME{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    length::Entry{T} = Entry(T, "m", "Some decription") # it's ok not to have a default, it forces users to think about what a parameter should be
    verbose::Entry{Bool} = Entry(Bool, "", "Some other decription"; default=true)
    switch::Switch{Symbol} = Switch(Symbol, [:option_a, :option_b], "", "user can only select one of these"; default=:option_a)
end

# Defintion of the actor structure
# NOTE: To be valid all actors must have `dd::IMAS.dd` and `par::FUSEparameters__ActorNAME`
mutable struct ActorNAME <: ???AbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorNAME  # Actors must carry with them the parameters they are run with
    something_else::??? # Some actors may want to carry something else with them
    # Inner constructor for the actor starting from `dd` and `par` (we generally refer to `par` as `act.ActorNAME`)
    # NOTE: Computation should not happen here since in workflows it is normal to instantiate
    #       an actor once `ActorNAME(dd, act.ActorNAME)` and then call `finalize(step(actor))` several times as needed.
    function ActorNAME(dd::IMAS.dd, par::FUSEparameters__ActorNAME; kw...)
        logging_actor_init(ActorNAME)
        par = par(kw...)
        return new(dd, par, something_else)
    end
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
   !!! tip Interactive documentation build
       One can call `include("make.jl")` over and over within the same Julia session to avoid dealing with startup time.

1. Check page by opening `FUSE/docs/build/index.html` page in web-browser.

1. To publish online, run in the `FUSE` folder:
   ```bash
   make web
   ```

!!! note
    Documentation files (PDF, DOC, XLS, PPT, ...) can be committed and pushed to the [FUSE_extra_files](https://github.com/ProjectTorreyPines/FUSE_extra_files) repository, and then linked directly from within the FUSE documentation, like this:

    ```markdown
    [video recording of the first FUSE tutorial](https://github.com/ProjectTorreyPines/FUSE_extra_files/raw/master/FUSE_tutorial_1_6Jul22.mp4)
    ```

## Examples

The `FUSE/examples` folder contains jupyter notebook that showcase some possible uses of FUSE.

!!! note
    When pushing changes to in a jupyter notebook, make sure that all the output cells are cleared 


## Tips and more

### Revise.jl

Install [Revise.jl](https://github.com/timholy/Revise.jl) to modify code and use the changes without restarting Julia.
We recommend adding `import Revise` to your `~/.julia/config/startup.jl` to automatically import Revise at the beginning of all Julia sessions.
All this can be done by running in the `FUSE` folder:

```bash
make revise
```

### Development in VSCode

[VSCode](https://code.visualstudio.com) is an excellent development environment for Julia.

FUSE uses the following VSCode settings for formatting the Julia code:

```json
{
    "files.autoSave": "onFocusChange",
    "workbench.tree.indent": 24,
    "editor.insertSpaces": true,
    "editor.tabSize": 4,
    "editor.detectIndentation": false,
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

!!! note
    To format Julia you will need to install `Julia Language Support` under the extensions tab (`<Command> + <shift> + x`)

### Tracking Julia precompilation

To see what is precompiled at runtime, you can add a Julia kernel with the `trace-compile` option to Jupyter

```julia
import IJulia
IJulia.installkernel("Julia tracecompile", "--trace-compile=stderr")
```

Then select the `Julia tracecompile` in jupyter-lab

## Running Julia within a Python environment

This can be particularly useful for benchmarking FUSE physics against existing Python routines (eg. in OMFIT)

1. Install `PyCall` in your Julia environment:

   ```
   export PYTHON="" # Sometimes one needs to empty the PYTHON environmental variable to install PyCall
   julia -e 'using Pkg; Pkg.add("PyCall"); Pkg.build("PyCall")'
   ```

   !!! note

       Python and Julia must be compiled for the same architecture.
       For example, to install Julia x64 in a Apple Silicon MACs:

       ```bash
       juliaup add release~x64
       export PYTHON=""
       julia +release~x64 -e 'using Pkg; Pkg.add("PyCall"); Pkg.build("PyCall")'
       ```

       You can make this verison your default one with

       ```bash
       juliaup default release~x64
       ```

1. Use pip to install the package PyJulia — remember to use the same Python passed to ENV["PYTHON"]:

   ```bash
   python3 –m pip install julia
   ```

1. Configure the communication between Julia and Python by running the following in the Python interpreter:

   ```python
   import julia
   julia.install()
   ```

   !!! note
       If you have more than one Julia version on our system, we could specify it with an argument:
       ```
       julia.install(julia="/Users/meneghini/.julia/juliaup/julia-1.8.5+0.x64.apple.darwin14/bin/julia")
       ```

1. Test the installation running the following in the Python interpreter run:

   ```python
   from julia import Main
   Main.eval('[x^2 for x in 0:4]')
   ```

1. Now, try something more useful:

   ```python
   from julia.api import Julia
   Julia(compiled_modules=False)
   def S(string): # from Python str to Julia Symbol
       return Main.eval(f"PyCall.pyjlwrap_new({string})")

   from julia import Main, IMAS, FUSE, Logging
   FUSE.logging(Logging.Info, actors=Logging.Debug);

   ini, act = FUSE.case_parameters(S(":FPP"), version=S(":v1_demount"), init_from=S(":scalars"), STEP=True);
   dd = FUSE.init(ini, act);

   eqt=dd.equilibrium.time_slice[-1]
   cp1d=dd.core_profiles.profiles_1d[-1]
   jFUSE = IMAS.Sauter_neo2021_bootstrap(eqt, cp1d, neo_2021=True)
   ```
