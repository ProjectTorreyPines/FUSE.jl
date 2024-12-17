if get(ENV, "FUSE_WITH_EXTENSIONS", "false") == "true"
    using Pkg
    Pkg.add("ThermalSystemModels")
    using ThermalSystemModels
end
using FUSE
using Test

@testset "warmup" begin
    println("== warmup ==")
    for round in (1, 2)
        dd = IMAS.dd()
        FUSE.warmup(dd)
    end
end

include("runtests_basics.jl")

include("runtests_cases.jl")

include("runtests_actors.jl")

include("runtests_workflow.jl")

include("runtests_init_expressions.jl")

if get(ENV, "FUSE_WITH_EXTENSIONS", "false") == "true"
    include("runtests_study.jl")
end

println(FUSE.timer)
