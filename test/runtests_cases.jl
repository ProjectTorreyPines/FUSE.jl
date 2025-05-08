using FUSE
using Test

@testset "test_cases" begin
    for testname in FUSE.available_test_cases()
        @testset "$(testname)" begin
            FUSE.TimerOutputs.reset_timer!(FUSE.timer, string(testname))
            FUSE.TimerOutputs.@timeit FUSE.timer string(testname) begin
                println("== $(testname) ==")
                dd = IMAS.dd()
                FUSE.test_case(testname, dd)
            end
        end
    end
end
