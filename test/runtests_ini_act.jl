using FUSE
using Test

@testset "ini_act" begin
    for (testname, (args, kw)) in FUSE.test_cases
        @testset "$testname" begin
            FUSE.TimerOutputs.reset_timer!(FUSE.timer, testname)
            FUSE.TimerOutputs.@timeit FUSE.timer "$testname" begin
                println("== $(testname) ==")
                FUSE.test_ini_act_save_load(args..., kw...)
            end
        end
    end
end
