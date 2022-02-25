function init_missing(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        for ids in keys(dd1)
            if length(keys(getproperty(dd, ids))) == 0
                setproperty!(dd, ids, getproperty(dd1, ids))
            end
        end
    end

    return dd
end