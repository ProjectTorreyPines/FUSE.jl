# Contributing

FUSE is a collaborative project that welcomes community contributions!

The `master` branch of ProjectTorreyPines repositories is write-protected. This means that even with write permissions to the repository, you'll not be able to push to `master` directly. Instead, we handle updates – be it new features or bug fixes – through branches and Pull Requests (PRs).

A crucial part of our PR process is **code review**. It is where your peers get to weigh in and ensure everything is up to standard before merging. When you create a PR, think about who on the team has the right expertise for the code you're working on, and assign them as reviewers. Their insights will not only help in maintaining code quality but also in catching any potential issues early. It is all about teamwork and making sure our code is the best it can be!

!!! note
    When working on a new feature that involves changes to FUSE and other ProjectTorreyPines repositories, you'll want to use the same branch name across these repositories. For example, if you're working on a branch named `my_new_feature` in both FUSE and IMAS, regression testing will be performed using the `my_new_feature` branches for FUSE and IMAS, along with the `master` branch of the other `ProjectTorreyPines` repositories.

## Add/modify entries in `dd`

The `dd` data structure is defined under the [IMASdd.jl](https://github.com/ProjectTorreyPines/IMASdd.jl) package. See the documentation there to how add/modify entries in `dd`.

## Write IMAS physics functions

IMAS physics and engineering functions are structured in [IMAS.jl](https://github.com/ProjectTorreyPines/IMAS.jl) under `IMAS/src/physics`. These functions use or modify the datastructure (dd) in some way and are used to calculate certain quantities or fill the data structure.

Let's say we want to create a function that calculates the DT fusion and then fill `core_sources` with the alpha heating from that source. Here is an example of writing it in a good way:

```julia
function DT_fusion_source!(dd::IMAS.dd)
    return DT_fusion_source!(dd.core_sources, dd.core_profiles)
end

"""
    DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)

Calculates DT fusion heating with an estimation of the alpha slowing down to the ions and electrons, modifies dd.core_sources
"""
function DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)
    # // actual implementation here //
end
```

!!! convention
    The documentation string is added to the specialized function `DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)` and the dispatch function `DT_fusion_source!(dd::IMAS.dd)` is added on top of the function

## Add/modify entries in `ini` and `act`

The functinoality of the `ini` and `act` parameters is implemented in the [SimulationParameters.jl](https://github.com/ProjectTorreyPines/SimulationParameters.jl) package.

* The `ini` parameters are all defined in the `FUSE/src/parameters_init.jl` file. Add/edit entries there.
* The `act` parameters of each actor are defined where that actor is defined. Add/edit entries there.

## Write a new actor

Actors are grouped into two main abstract types:

```julia
abstract type CompoundAbstractActor{D,P} <: AbstractActor{D,P} end
abstract type SingleAbstractActor{D,P} <: AbstractActor{D,P} end
```

CompoundAbstractActors are for actors that compound multiple actors underneath and are initalized with ```ActorNAME(dd, par, act)``` while SingleAbstractActors are single actors initalized with ```ActorNAME(dd, par)```

The definition of each FUSE actor follows a well defined pattern.
**DO NOT** deviate from this pattern. This is important to ensure modularity and compostability of the actors.

```julia
# Definition of the `act` parameters relevant to the actor
Base.@kwdef mutable struct FUSEparameters__ActorNAME{T<:Real} <: ParametersActor{T}
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

## Add a new material 

Material properties for supported fusion-relevant materials are stored in the [FusionMaterials.jl](https://github.com/ProjectTorreyPines/FusionMaterials.jl) package, specifically in `FusionMaterials/src/materials.jl`. Properties of each material can be accessed by calling the `Material` function with the material name as a symbol passed as the function argument. 

To add a new material whose properties can be accessed in FUSE, first add a function to `materials.jl` called Material with the function argument being your material's name. In the body of the function, assign the material's name (as a string, all lowercase, and with any spaces filled by underscores), type (as a list containing each possible IMAS BuildLayerType the material could be assigned to), density (in `kg/m^3`) and unit cost (in US dollars per kilogram). Include a comment providing a link to the source from which the unit cost was taken. 

Below is an example of a complete Material function for a non-superconductor material (more about superconductor materials below): 

```julia 
function Material(::Type{Val{:graphite}};)
	mat = Material()
	mat.name = "graphite" # string with no spaces
	mat.type = [IMAS._wall_] # list of allowable layer types for this material
	mat.density = 1.7e3 # in kg/m^3
	mat.unit_cost = 1.3 # in US$/kg, include source as a comment # source: https://businessanalytiq.com/procurementanalytics/index/graphite-price-index/
	return mat
end
```

If the material is a superconductor that is meant to be assigned to magnet-type layers, additional characteristics need to be defined. First, add the relevant critical current density scaling for the chosen superconductor material as a function in `FusionMaterials/src/jcrit.jl`. Then, assign the technology parameters for the material (temperature, steel fraction, void fraction, and ratio of superconductor to copper) to their respective fields in coil_tech within the coil_technology function in `FUSE/src/technology.jl`. Finally, call the critical current density scaling function within the newly written Material function in `materials.jl` and assign the output critical current density and critical magnetic field to the material object. The coil_tech object should be passed as an argument to the Material function, along with the external B field, and used to calculate the critical current density and critical magnetic field. 

Below is an example of a complete superconductor Material function: 

```julia 
function Material(::Type{Val{:rebco}}; coil_tech::Union{Missing, IMAS.build__pf_active__technology, IMAS.build__oh__technology, IMAS.build__tf__technology} = missing, Bext::Union{Real, Missing} = missing)
	mat = Material()
	mat.name = "rebco"
	mat.type = [IMAS._tf_, IMAS._oh_]
	mat.density = 6.3
	mat.unit_cost = 7000

	if !ismissing(coil_tech)
		Jcrit_SC, Bext_Bcrit_ratio = ReBCO_Jcrit(Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature) # A/m^2
		fc = fraction_conductor(coil_tech)
		mat.critical_current_density = Jcrit_SC * fc
		mat.critical_magnetic_field = Bext / Bext_Bcrit_ratio
	end

	return mat
end
```

The function `ReBCO_Jcrit` is the critical current density function for this material. 

You can then access the parameters of your material by calling the function you've created. For example, access the material's density anywhere in FUSE by calling: 

```julia
my_mat_density = Material(:my_mat).density
```

## Profiling and writing fast Julia code

First let's do some profiling to identify problemetic functions:

1. Run FUSE from the VScode Julia REPL (`<Command-Shift-p>` and then `Julia: Start REPL`)
   ```julia
   using FUSE
   FUSE.logging(Logging.Info; actors=Logging.Info);
   ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars, STEP=true);
   dd = IMAS.dd()
   FUSE.init(dd, ini, act)
   FUSE.ActorWholeFacility(dd, act); # call this once to precompile
   ```

   !!! note
   Alternatively one can create a `profile.jl` file in the `FUSE/playground` folder, write Julia code in that file, select the code to execute and run it with `<Control-Return>`.

1. Use `@time` to monitor execution time and allocations

1. For functions that return very quickly one can use `BenchmarkTooks.@benchmark`

1. Graphical profiling of the execution time is a powerful way to understand where time is spent
   ```julia
   @profview FUSE.ActorWholeFacility(dd, act);
   ```
   where `FUSE.ActorWholeFacility(dd, act);` can really be any function that we care about

1. Look at allocations
   ```julia
   @profview_allocs FUSE.ActorWholeFacility(dd, act);
   ```

1. We can decide how finely to comb for bottlenecks by setting `sample_rate` in `@profview` and `@profview_allocs`:
   ```
   @profview_allocs f(args...) [sample_rate=0.0001] [C=false]
   ```

!!! note
    To move forward we have to [understand how to write performant Julia code](https://docs.julialang.org/en/v1/manual/performance-tips).

Let's now investigate where the issue is with the function that we have identified. For this we have several tools at our disposal:

* [`@code_warntype`](https://docs.julialang.org/en/v1/manual/performance-tips/#man-code-warntype): static analyzer built-in with Julia
  * only looks at types that are inferred at runtime
  * reports types only for the target function
  * `@code_warntype function()`

* [JET](https://github.com/aviatesk/JET.jl): static analyzer integrated with VScode
  * can detect different possible issues, including types inferred at runtime
  * JET goes deep into functions
  * `JET.@report_opt function()` reports dynamic dispatch
  * `JET.@report_call function()` reports type errors
  * `JET.@report_call target_modules=(FUSE,IMAS,IMAS.IMASdd, ) FUSE.ActorNeutronics(dd,act);`

* [Cthulhu](https://github.com/JuliaDebug/Cthulhu.jl): interactive static analyzer
  * `Cthulhu.@descend function()`

## Build the documentation

1. To build the documentation, in the `FUSE/docs` folder, start Julia then:
   ```julia
   ] activate .
   include("make.jl")
   ```
   !!! tip Interactive documentation build
       One can call `include("make.jl")` over and over within the same Julia session to avoid dealing with startup time.

1. Check page by opening `FUSE/docs/build/index.html` page in web-browser.

1. The online documentation is built after each commit to `master` via GitHub actions.

!!! note
    Documentation files (PDF, DOC, XLS, PPT, ...) can be committed and pushed to the [FUSE\_extra\_files](https://github.com/ProjectTorreyPines/FUSE_extra_files) repository, and then linked directly from within the FUSE documentation, like this:

    ```markdown
    [video recording of the first FUSE tutorial](https://github.com/ProjectTorreyPines/FUSE_extra_files/raw/master/FUSE_tutorial_1_6Jul22.mp4)
    ```

## Add/modify examples

The [FuseExamples repository](https://github.com/ProjectTorreyPines/FuseExamples) contains jupyter notebook that showcase some possible uses of FUSE.

!!! note
    When committing changes to in a jupyter notebook, make sure that all the output cells are cleared! This is important to keep the size of the repository in check.

## Use Revise.jl

Install [Revise.jl](https://github.com/timholy/Revise.jl) to modify code and use the changes without restarting Julia.
We recommend adding `import Revise` to your `~/.julia/config/startup.jl` to automatically import Revise at the beginning of all Julia sessions. This can be done by running:

```bash
fusebot install_revise
```

## Develop in VScode

[VScode](https://code.visualstudio.com) is an excellent development environment for Julia.

FUSE uses the following VScode settings for formatting the Julia code:

```json
{
    "files.autoSave": "onFocusChange",
    "workbench.tree.indent": 24,
    "editor.insertSpaces": true,
    "editor.tabSize": 4,
    "editor.detectIndentation": false,
    "[julia]": {
        "editor.defaultFormatter": "julialang.language-julia"
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
    To add these settings to VScode add these lines to: `<Command> + <shift> + p` -> `Preferences: Open User Settings (JSON)`

!!! note
    To format Julia you will need to install `Julia Language Support` under the extensions tab (`<Command> + <shift> + x`)

## Track Julia precompilation

To see what is precompiled at runtime, you can add a Julia kernel with the `trace-compile` option to Jupyter

```julia
import IJulia
IJulia.installkernel("Julia tracecompile", "--trace-compile=stderr")
```

Then select the `Julia tracecompile` in jupyter-lab


!!! note
    If you want to remove jupyter kernels you don't use anymore you can list them first with ```jupyter kernelspec list``` and remove via ```jupyter kernelspec remove <your kernel>```

## Run Julia within a Python environment

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
