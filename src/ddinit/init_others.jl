function init_missing(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        for ids in keys(dd1)
            if length(keys(getproperty(dd, ids))) == 0
                data = getproperty(dd1, ids)
                dd.global_time = max(dd.global_time, maximum(data.time))
                setproperty!(dd, ids, data)
            end
        end
    end

    return dd
end