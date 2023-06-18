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
            ini, act = case_parameters(:FPP; version=:v1_demount, init_from=:scalars, STEP=:true)
            init(dd, ini, act)
        end
        act.ActorEquilibriumTransport.max_iter = 1
        ActorWholeFacility(dd, act)
        return IMAS.freeze(dd)
    end
end
