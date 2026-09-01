using FUSE
using Test

# By default only run a subset of the test cases, chosen to cover the different
# code paths (init from ods/scalars/experiment/database, H/L-mode, time-dependent)
# rather than every machine. KDEMO is covered by runtests_warmup.jl.
# The nightly scheduled CI run sets FUSE_ALL_TEST_CASES=true and runs everything.
const default_test_cases = [:ITER_ods, :ITER_scalars, :ITER_time, :D3D_Hmode, :D3D_Lmode, :JET_HDB5, :FPP]

@testset "test_cases" begin
    if get(ENV, "FUSE_ALL_TEST_CASES", "false") == "true"
        test_cases = FUSE.available_test_cases()
    else
        test_cases = default_test_cases
    end
    for testname in test_cases
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
