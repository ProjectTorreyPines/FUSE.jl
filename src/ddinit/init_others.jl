"""
    init_missing_from_ods(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize missing IDSs from ODS, only if `ini.general.init_from == :ods`.
"""
function init_missing_from_ods(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
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
    end

    # target
    for field in [:power_electric_net, :flattop_duration, :tritium_breeding_ratio, :cost]
        value = getproperty(ini.target, field, missing)
        if value !== missing
            setproperty!(dd.target, field, value)
        end
    end

    return dd
end