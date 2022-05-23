"""
    init_missing_from_ods(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)

Initialize missing IDSs from ODS, only if `ini.general.init_from == :ods`
"""
function init_missing_from_ods(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)
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

    return dd
end