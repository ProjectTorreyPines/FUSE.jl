"""
    init_currents!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.core_profiles` and `dd.core_sources` ohmic and bootstrap currents and sources starting from `ini` and `act` parameters
"""
function init_currents!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_currents")
    TimerOutputs.@timeit timer "init_currents" begin
        init_from = ini.general.init_from

        if (init_from == :scalars) || ismissing(dd.core_profiles.profiles_1d[], :j_ohmic)
            ActorCurrent(dd, act; ip_from=:pulse_schedule, vloop_from=:pulse_schedule)
            IMAS.ohmic_source!(dd)
        end

        return dd
    end
end
