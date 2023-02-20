# ******************************************
# save/load simulation
# ******************************************
"""
    save(
        dd::IMAS.dd,
        ini::ParametersAllInits,
        act::ParametersAllActors,
        dirname::AbstractString;
        freeze::Bool=true,
        format::Symbol=:hdf)

Save FUSE dd, ini, act files in a folder
`dd` can be saved in JSON or HDF format
"""
function save(
    dd::IMAS.dd,
    ini::ParametersAllInits,
    act::ParametersAllActors,
    dirname::AbstractString;
    freeze::Bool=true,
    format::Symbol=:hdf)

    @assert format in [:hdf, :json] "format must be either `:hdf` or `:json`"
    mkdir(dirname) # purposely error if directory exists or path does not exist
    if format == :hdf
        IMAS.imas2hdf(dd, joinpath(dirname, "dd.h5"); freeze)
    elseif format == :json
        IMAS.imas2json(dd, joinpath(dirname, "dd.json"); freeze)
    end
    ini2json(ini, joinpath(dirname, "ini.json"))
    act2json(act, joinpath(dirname, "act.json"))
    return nothing
end

"""
    load(dirname::AbstractString)

Read dd, ini, act from files in a folder
"""
function load(dirname::AbstractString)
    if isfile(joinpath(dirname, "dd.h5"))
        dd = IMAS.hdf2imas(joinpath(dirname, "dd.h5"))
    else
        dd = IMAS.json2imas(joinpath(dirname, "dd.json"))
    end
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
    return load(dir, [extract])[1]
end

function load(dir::AbstractString, extracts::Vector{Dict{Symbol,Function}})::Vector{Dict{Symbol,Any}}
    dd, ini, act = FUSE.load(dir)
    out = Dict{Symbol,Any}[]
    for extract in extracts
        results = Dict{Symbol,Any}()
        for key in collect(keys(extract))
            value = try
                extract[key](dd, ini, act)
            catch e
                NaN
            end
            results[key] = value
        end
        push!(out, results)
    end
    return out
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
    return load(dirs, [extract])[1]
end

function load(dirs::AbstractVector{<:AbstractString}, extracts::Vector{Dict{Symbol,Function}})::Vector{DataFrames.DataFrame}
    dfs = DataFrames.DataFrame[]
    for kdf in 1:length(extracts)
        df = DataFrames.DataFrame(load(dirs[1], extracts[kdf]))
        for k in 2:length(dirs)
            push!(df, df[1, :])
        end
        push!(dfs, df)
    end
    p = ProgressMeter.Progress(length(dirs); showspeed=true)
    GC.enable(false) #https://github.com/JuliaIO/HDF5.jl/pull/1049
    try
        Threads.@threads for k in 1:length(dirs)
            tmp = load(dirs[k], extracts)
            for kdf in 1:length(extracts)
                dfs[kdf][k, :] = tmp[kdf]
            end
            ProgressMeter.next!(p)
        end
    finally
        GC.enable(true)
    end
    dfs
end
