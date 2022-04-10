function init_missing(dd::IMAS.dd, ini::InitParameters, act::ActorParameters)
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