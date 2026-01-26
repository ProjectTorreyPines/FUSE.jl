# Contributing

FUSE is a collaborative project, and we welcome contributions from the community!

## Setting up a developer environment

Follow these steps to set up your development environment for contributing to the FUSE codebase (or any other Julia package).

1. Install and configure Revise.jl

   Install [Revise.jl](https://github.com/timholy/Revise.jl) to modify code and see changes without restarting Julia. We recommend adding `import Revise` to your `~/.julia/config/startup.jl` to automatically load Revise in all Julia sessions. Run the following command:

   ```bash
   fusebot install_revise
   ```

   Restart your Julia sessions and Jupyter kernels to apply this change.

2. By default, development packages are located in `~/.julia/dev`. Since accessing hidden folders like `~/.julia` can be cumbersome, you may find it helpful to create a symbolic link to this dev folder (you can name the link `~/julia_dev` or something more convenient):

   ```bash
   ln -s ~/.julia/dev ~/julia_dev
   ```

3. Instead of manually cloning repositories, we recommend using Julia's `Pkg` system. Here’s a quick guide:

   ```julia
   ] dev FUSE # clones the `master` branch of FUSE
   ] dev FUSE\#feature_branch # clones a specific feature branch
   ] dev FUSE@v1.2.3 # checks out a specific release version
   ] free FUSE # stops development mode, reverting to standard usage
   ] status # shows the status of all packages in use and development
   ] up # upgrades all packages (note: manage `dev` packages manually)
   ```

   These commands work not only for FUSE but for any Julia package, including [packages that FUSE depends on](https://fuse.help/dev/deps.html).

    Familiarize yourself with [managing Julia packages](https://pkgdocs.julialang.org/v1/managing-packages/), as this will be helpful when contributing to FUSE or other Julia projects.

3. Once set up, you can open `~/julia_dev` in VS Code, make changes to the code, and immediately see the effects of your changes live in your Julia sessions and Jupyter notebooks (thanks to `Revise.jl`).

## Git branches

The `master` and `dev` branches of ProjectTorreyPines repositories are write-protected. This means that even with write permissions to the repository, you'll not be able to push to them directly. Instead, we handle updates – be it new features or bug fixes – through branches and Pull Requests (PRs).

A crucial part of our PR process is **code review**. It is where your peers get to weigh in and ensure everything is up to standard before merging. When you create a PR, think about who on the team has the right expertise for the code you're working on, and assign them as reviewers. Their insights will not only help in maintaining code quality but also in catching any potential issues early. It is all about teamwork and making sure our code is the best it can be!

!!! note
    When working on a new feature that involves changes to FUSE and other ProjectTorreyPines repositories, you'll want to use the same branch name across these repositories. For example, if you're working on a branch named `my_new_feature` in both FUSE and IMAS, regression testing will be performed using the `my_new_feature` branches for FUSE and IMAS, along with the `master` branch of the other `ProjectTorreyPines` repositories.

## Add/modify entries in `dd`

The `dd` data structure is defined under the [IMASdd.jl](https://github.com/ProjectTorreyPines/IMASdd.jl) package.
See the documentation there for how to add/modify entries in `dd`.

## Write a new actor

See the Jupyter [new actor tutorial](https://github.com/ProjectTorreyPines/FuseExamples/blob/master/new_actor.ipynb)

## Add/modify examples

The [FuseExamples repository](https://github.com/ProjectTorreyPines/FuseExamples) contains jupyter notebooks that showcase some possible uses of FUSE.

!!! note
    When committing changes to a jupyter notebook, make sure that all the output cells are cleared! This is important to keep the size of the repository in check.

## Add/modify a new cluster

The `parallel_environment()` function is used to easily setup the working environment for distributed jobs. Edit it [here](https://github.com/search?q=repo%3AProjectTorreyPines%2FFUSE.jl+%22function+parallel_environment%22&type=code) with your own cluster, following the others as a template.

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
    The documentation string is added to the specialized function `DT_fusion_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)` and the dispatch function `DT_fusion_source!(dd::IMAS.dd)` is added on above

## Add/modify entries in `ini` and `act`

The functionality of the `ini` and `act` parameters is implemented in the [SimulationParameters.jl](https://github.com/ProjectTorreyPines/SimulationParameters.jl) package.

* The `ini` parameters are all defined in the `FUSE/src/parameters_init.jl` file. Add/edit entries there.
* The `act` parameters of each actor are defined where that actor is defined. Add/edit entries there.

## Profiling and writing fast Julia code

First let's do some profiling to identify problematic functions:

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

1. For functions that return very quickly one can use `BenchmarkTools.@benchmark`

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
       If you have more than one Julia version on your system, we could specify it with an argument:
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
