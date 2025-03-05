"""
    init_missing_from_ods!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize missing IDSs from ODS, only if `ini.general.init_from == :ods`.
"""
function init_missing_from_ods!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_missing_from_ods")
    TimerOutputs.@timeit timer "init_missing_from_ods" begin
        init_from = ini.general.init_from

        if init_from == :ods
            for field in keys(dd1)
                ids1 = getproperty(dd1, field)
                ids = getproperty(dd, field)
                if !IMAS.hasdata(ids) && IMAS.hasdata(ids1)
                    TimerOutputs.reset_timer!("init_$(field)_from_ods")
                    TimerOutputs.@timeit timer "init_$(field)_from_ods" begin
                        setproperty!(dd, field, deepcopy(ids1))
                    end
                end
            end
        end

        return dd
    end
end

"""
    init_requirements!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize dd.requirements `ini.requirements`
"""
function init_requirements!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_requirements")
    TimerOutputs.@timeit timer "init_requirements" begin

        # NOTE: `log10_flattop_duration` (used when running optimizations) wins over `flattop_duration`
        for field in sort([field for field in fieldnames(typeof(ini.requirements)) if string(field)[1] != '_']; by=x -> startswith(string(x), "log10_"))
            value = getproperty(ini.requirements, field, missing)
            if value !== missing
                if startswith(string(field), "log10_")
                    field = Symbol(replace(string(field), "log10_" => ""))
                    value = 10.0 .^ value
                end
                setproperty!(dd.requirements, field, value)
            end
        end

        return dd
    end
end
