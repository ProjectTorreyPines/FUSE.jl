"""
    init_missing_from_ods(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize missing IDSs from ODS, only if `ini.general.init_from == :ods`.
"""
function init_missing_from_ods(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    TimerOutputs.reset_timer!("init_missing_from_ods")
    TimerOutputs.@timeit timer "init_missing_from_ods" begin
        init_from = ini.general.init_from

        if init_from == :ods
            dd1 = IMAS.json2imas(ini.ods.filename)
            for ids in keys(dd1)
                if length(keys(getproperty(dd, ids))) == 0
                    data = getproperty(dd1, ids)
                    if !ismissing(data, :time)
                        dd.global_time = max(dd.global_time, maximum(data.time))
                    end
                    setproperty!(dd, ids, data)
                end
            end
            if :core_profiles ∈ keys(dd1)
                ne_ped, w_ped = IMAS.pedestal_finder(dd.core_profiles.profiles_1d[].electrons.density_thermal, dd.core_profiles.profiles_1d[].grid.psi_norm)
                ped_summ = dd.summary.local.pedestal
                cp1d = dd.core_profiles.profiles_1d[]
                @ddtime ped_summ.n_e.value = ne_ped
                @ddtime ped_summ.position.rho_tor_norm = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.grid.rho_tor_norm).(1 - w_ped)
                @ddtime ped_summ.zeff.value = IMAS.interp1d(dd.core_profiles.profiles_1d[].grid.rho_tor_norm, dd.core_profiles.profiles_1d[].zeff).(1 - w_ped)
            end
        end

        # requirements
        # NOTE: `log10_flattop_duration` (used when running optimizations) wins over `flattop_duration`
        for field in sort([string(field) for field in fieldnames(typeof(ini.requirements)) if string(field)[1]!='_'], by=x->startswith(x,"log10_"))
            value = getproperty(ini.requirements, field, missing)
            if value !== missing
                if startswith(string(field), "log10_")
                    field = Symbol(replace(string(field), "log10_" => ""))
                    value = 10.0.^value
                end
                setproperty!(dd.requirements, field, value)
            end
        end

        return dd
    end
end
