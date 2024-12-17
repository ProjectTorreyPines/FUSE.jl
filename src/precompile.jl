"""
    warmup(dd::IMAS.dd)

Function used to precompile the majority of FUSE

NOTE: This function is also useful for profiling
"""
function warmup(dd::IMAS.dd)
    empty!(dd)
    TimerOutputs.reset_timer!("warmup")
    return TimerOutputs.@timeit timer "warmup" begin
        ini, act = case_parameters(:KDEMO)
        ini_act_tests_customizations!(ini, act)

        TimerOutputs.@timeit timer "init" begin
            init(dd, ini, act)
        end

        TimerOutputs.@timeit timer "ActorWholeFacility" begin
            ActorWholeFacility(dd, act)
        end

        TimerOutputs.@timeit timer "freeze" begin
            frozen_dd = IMAS.freeze(dd)
        end

        return frozen_dd
    end
end
