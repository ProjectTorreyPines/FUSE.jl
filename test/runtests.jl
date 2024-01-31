using FUSE
using Test

@testset "warmup" begin
    for round in (1, 2)
        dd = IMAS.dd()
        FUSE.warmup(dd)
    end
end

include("runtests_basics.jl")

include("runtests_cases.jl")

include("runtests_actors.jl")

include("runtests_workflow.jl")

println(FUSE.timer)
