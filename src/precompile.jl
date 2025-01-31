"""
    warmup(dd::IMAS.dd)

Function used to precompile the majority of FUSE
"""
function warmup(dd::IMAS.dd)
    empty!(dd)

    ini, act = case_parameters(:KDEMO)

    ini_act_tests_customizations!(ini, act)

    init(dd, ini, act)

    ActorWholeFacility(dd, act)

    TimerOutputs.@timeit timer "freeze" begin
        frozen_dd = IMAS.freeze(dd)
    end

    return frozen_dd
end
