"""
    warmup(dd::IMAS.dd)

Function used to precompile the majority of FUSE

This function is also useful for profiling
"""
function warmup(dd::IMAS.dd)
    empty!(dd)
    TimerOutputs.reset_timer!("warmup")
    return TimerOutputs.@timeit timer "warmup" begin
        TimerOutputs.@timeit timer "init" begin
            ini, act = case_parameters(:FPP)
            init(dd, ini, act)
        end
        act.ActorStationaryPlasma.max_iter = 1
        ActorWholeFacility(dd, act)
        return IMAS.freeze(dd)
    end
end
