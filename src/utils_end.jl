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
    return dirname
end

"""
    load(dirname::AbstractString)

Returns (dd, ini, act) from files read in a folder

`dd` can be in in JSON `dd.json` or HDF `dd.h5` format.

Returns `missing` for files are not there
"""
function load(dirname::AbstractString)
    if isfile(joinpath(dirname, "dd.h5"))
        dd = IMAS.hdf2imas(joinpath(dirname, "dd.h5"))
    elseif isfile(joinpath(dirname, "dd.json"))
        dd = IMAS.json2imas(joinpath(dirname, "dd.json"))
    else
        dd = missing
    end
    if isfile(joinpath(dirname, "ini.json"))
        ini = json2ini(joinpath(dirname, "ini.json"))
    else
        ini = missing
    end
    if isfile(joinpath(dirname, "act.json"))
        act = json2act(joinpath(dirname, "act.json"))
    else
        act = missing
    end
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

function load(dirs::AbstractVector{<:AbstractString}, extracts::Vector{Dict{Symbol,Function}}; filter_invalid::Bool=true)::Vector{DataFrames.DataFrame}
    # at first we load the data all in the same dataframe (for filtering)
    all_extracts = Dict{Symbol,Function}()
    for extract in extracts
        merge!(all_extracts, extract)
    end

    # allocate memory
    df = DataFrames.DataFrame(load(dirs[1], all_extracts))
    for k in 2:length(dirs)
        push!(df, df[1, :])
    end
    
    # load the data
    p = ProgressMeter.Progress(length(dirs); showspeed=true)
    GC.enable(false) #https://github.com/JuliaIO/HDF5.jl/pull/1049
    try
        Threads.@threads for k in 1:length(dirs)
            tmp = load(dirs[k], all_extracts)
            df[k, :] = tmp
            ProgressMeter.next!(p)
        end
    finally
        GC.enable(true)
    end

    # filter
    if filter_invalid
        df = filter(row -> all(x -> !(x isa Number && isnan(x)), row), df)
    end

    # split into separate dataframes
    dfs = DataFrames.DataFrame[]
    for extract in extracts
        push!(dfs, df[:,collect(keys(extract))])
    end

    return dfs
end
