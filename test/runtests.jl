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
    end

    include("runtests_warmup.jl")

    include("runtests_basics.jl")

    include("runtests_cases.jl")

    include("runtests_actors.jl")

    include("runtests_init_expressions.jl")

    # Conditionally run extension tests if environment variable is set.
    if get(ENV, "FUSE_WITH_EXTENSIONS", "false") == "true"
        include("runtests_study.jl")
    end

    println(FUSE.timer)
end
