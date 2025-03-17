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

function warmup()

    println("Run ActorWholeFacility for FPP")
    dd = IMAS.dd()
    actor_logging(dd, false)
    ini, act = case_parameters(:FPP)
    ini_act_tests_customizations!(ini, act)
    init(dd, ini, act)
    ActorWholeFacility(dd, act)
    IMAS.freeze(dd)
    println()

    println("Run ActorStationaryPlasma for ITER, :ods")
    empty!(dd)
    actor_logging(dd, false)
    ini, act = case_parameters(:ITER, init_from=:ods)
    ini_act_tests_customizations!(ini, act)
    act.ActorEquilibrium.model = :FRESCO
    act.ActorCurrent.model = :SteadyStateCurrent;
    init(dd, ini, act)
    ActorStationaryPlasma(dd, act)
    println()

    println("Run ActorDynamicPlasma for D3D, :default")
    empty!(dd)
    actor_logging(dd, false)
    ini, act = case_parameters(:D3D, :default)
    ini_act_tests_customizations!(ini, act)
    act.ActorEquilibrium.model = :FRESCO
    act.ActorDynamicPlasma.Δt = 0.01
    act.ActorDynamicPlasma.Nt = 2
    act.ActorDynamicPlasma.ip_controller = true
    act.ActorDynamicPlasma.time_derivatives_sources = true
    init(dd, ini, act)
    ActorDynamicPlasma(dd, act)

    return IMAS.freeze(dd)
end
