"""
    warmup()

Function used to precompile the majority of FUSE
"""
function warmup()
    dd = IMAS.dd()
    return warmup(dd)
end

function warmup(dd::IMAS.dd)
    TimerOutputs.reset_timer!("warmup")
    return TimerOutputs.@timeit timer "warmup" begin
        TimerOutputs.@timeit timer "init" begin
            ini, act = case_parameters(:FPP; version=:v1_demount, init_from=:scalars, STEP=:true)
            init(dd, ini, act)
        end
        ActorWholeFacility(dd, act)
        return IMAS.freeze(dd)
    end
end
