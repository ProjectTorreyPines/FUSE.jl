"""
    init_currents(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.core_profiles` and `dd.core_sources` ohmic and bootstrap currents and sources starting from `ini` and `act` parameters
"""
function init_currents(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    TimerOutputs.reset_timer!("init_currents")
    TimerOutputs.@timeit timer "init_currents" begin
        init_from = ini.general.init_from
        dd.global_time = ini.time.simulation_start

        if (init_from == :scalars) || ismissing(dd.core_profiles.profiles_1d[], :j_ohmic)
            ActorSteadyStateCurrent(dd, act)
        end

        return dd
    end
end
