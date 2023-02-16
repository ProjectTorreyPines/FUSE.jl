# ******************************************
# save/load simulation
# ******************************************
"""
    save(
        dd::IMAS.dd,
        ini::ParametersAllInits,
        act::ParametersAllActors,
        dirname::AbstractString)

Save FUSE dd, ini, act as JSON files in a folder
"""
function save(
    dd::IMAS.dd,
    ini::ParametersAllInits,
    act::ParametersAllActors,
    dirname::AbstractString)

    mkdir(dirname)
    IMAS.imas2json(dd, joinpath(dirname,"dd.json"))
    ini2json(ini,joinpath(dirname,"ini.json"))
    act2json(act,joinpath(dirname,"act.json"))
    return nothing
end
