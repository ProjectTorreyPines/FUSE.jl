# Check if specific test files are requested via ARGS
if !isempty(ARGS)
    for testfile in ARGS
        @info "Running test file: $testfile"
        include(testfile)
    end
else
    # Default behavior: run all tests

    # Conditionally run extension tests if environment variable is set.
    if get(ENV, "FUSE_WITH_EXTENSIONS", "false") == "true"
        using Pkg
        Pkg.add("ThermalSystemModels")
        using ThermalSystemModels
        Pkg.add(url="https://github.com/OpenFUSIONToolkit/GPEC.git", rev="bfc60610ee04511c830a502c704f269671714691")
        using GeneralizedPerturbedEquilibrium
    end

    include("runtests_warmup.jl")

    include("runtests_basics.jl")

    include("runtests_simplenb_energy_components.jl")
    include("runtests_sawteeth_flatten_bootstrap.jl")
    include("runtests_qed_sawtooth_crash.jl")
    include("runtests_simplenb_axis_regularization.jl")
    include("runtests_simplenb_power_fractions.jl")

    include("runtests_cases.jl")

    include("runtests_actors.jl")

    include("runtests_init_expressions.jl")

    # Conditionally run extension tests if environment variable is set.
    if get(ENV, "FUSE_WITH_EXTENSIONS", "false") == "true"
        include("runtests_study.jl")
        include("runtests_gpec.jl")
    end

    include("runtests_study_database.jl")

    # forward-mode AD: whole actors on a dual-typed dd, every derivative against finite differences
    include("runtests_dual_actors.jl")
    
    println(FUSE.timer)
end
