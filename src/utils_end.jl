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
    dirname::AbstractString;
    freeze::Bool=true)

    mkdir(dirname) # will purposly error if directory exists or path to directory does not exist
    IMAS.imas2json(dd, joinpath(dirname, "dd.json"); freeze)
    ini2json(ini, joinpath(dirname, "ini.json"))
    act2json(act, joinpath(dirname, "act.json"))
    return nothing
end

"""
    load(dirname::AbstractString)

Read dd, ini, act from JSON files in a folder
"""
function load(dirname::AbstractString)
    dd = IMAS.json2imas(joinpath(dirname, "dd.json"))
    ini = json2ini(joinpath(dirname, "ini.json"))
    act = json2act(joinpath(dirname, "act.json"))
    return dd, ini, act
end

"""
    load(dir::AbstractString, extract::Dict{Symbol,Function})::Dict{Symbol,Any}

Read dd, ini, act from JSON files in a folder and extract some data from them

`extract` functions should accept `dd, ini, act` as inputs, like this:

    extract = Dict(
            :beta_normal => (dd,ini,act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
            :time => (dd,ini,act) -> @ddtime(dd.equilibrium.time),
            :R0 => (dd,ini,act) -> ini.equilibrium.R0
        )

"""
function load(dir::AbstractString, extract::Dict{Symbol,Function})::Dict{Symbol,Any}
    results = Dict{Symbol,Any}()
    dd, ini, act = FUSE.load(dir)
    for key in collect(keys(extract))
        results[key] = extract[key](dd, ini, act)
    end
    results
end

"""
    load(dirs::AbstractVector{<:AbstractString}, extract::Dict{Symbol,Function})::DataFrames.DataFrame

Read dd, ini, act from JSON files from multiple directores and extract some data from them

`extract` functions should accept `dd, ini, act` as inputs, like this:

    extract = Dict(
            :beta_normal => (dd,ini,act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
            :time => (dd,ini,act) -> @ddtime(dd.equilibrium.time),
            :R0 => (dd,ini,act) -> ini.equilibrium.R0
        )
"""
function load(dirs::AbstractVector{<:AbstractString}, extract::Dict{Symbol,Function})::DataFrames.DataFrame
    df = DataFrames.DataFrame(load(dirs[1], extract))
    for k in 1:length(dirs)-1
        push!(df, df[1, :])
    end
    p = ProgressMeter.Progress(length(dirs); showspeed=true)
    Threads.@threads for (k, dir) in collect(enumerate(dirs))
        df[k, :] = load(dirs[k], extract)
        ProgressMeter.next!(p)
    end
    df
end
