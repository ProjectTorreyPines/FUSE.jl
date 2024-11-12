using FUSE
using Test

@testset "use_cases" begin
    for (testname, (args, kw)) in FUSE.test_cases
        @testset "$testname" begin
            FUSE.TimerOutputs.reset_timer!(FUSE.timer, testname)
            FUSE.TimerOutputs.@timeit FUSE.timer "$testname" begin
                println("== $(testname) ==")
                dd = IMAS.dd()
                FUSE.test(dd, args...; kw...)
            end
        end
    end
end
