"""
    warmup()

Function used to precompile the majority of FUSE
"""
function warmup()
    dd = IMAS.dd()
    return warmup(dd)
end

function warmup(dd::IMAS.dd)
    TimerOutputs.reset_timer!(to, "warmup")
    return TimerOutputs.@timeit to "warmup" begin
        TimerOutputs.@timeit to "init" begin
            ini, act = case_parameters(:FPP; version=:v1_demount, init_from=:scalars)
            init(dd, ini, act)
        end
        ActorWholeFacility(dd, act)
        return IMAS.freeze(dd)
    end
end

if VERSION > v"1.8"
    SnoopPrecompile.@precompile_setup begin
        SnoopPrecompile.@precompile_all_calls begin
            FUSE.warmup()
        end
    end
end