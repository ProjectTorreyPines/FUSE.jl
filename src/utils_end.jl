# ******************************************
# save/load simulation
# ******************************************
"""
    save(
        dd::IMAS.dd,
        ini::ParametersAllInits,
        act::ParametersAllActors,
        savedir::AbstractString;
        freeze::Bool=true,
        format::Symbol=:json)

Save FUSE dd, ini, act files in a folder

`dd` can be saved in JSON or HDF format
"""
function save(
    savedir::AbstractString,
    dd::IMAS.dd,
    ini::ParametersAllInits,
    act::ParametersAllActors;
    freeze::Bool=true,
    format::Symbol=:json)

    @assert format in [:hdf, :json] "format must be either `:hdf` or `:json`"
    mkdir(savedir) # purposely error if directory exists or path does not exist
    if format == :hdf
        IMAS.imas2hdf(dd, joinpath(savedir, "dd.h5"); freeze)
    elseif format == :json
        IMAS.imas2json(dd, joinpath(savedir, "dd.json"); freeze)
    end
    ini2json(ini, joinpath(savedir, "ini.json"))
    act2json(act, joinpath(savedir, "act.json"))
    return savedir
end

"""
    save(
        savedir::AbstractString,
        dd::IMAS.dd,
        ini::ParametersAllInits,
        act::ParametersAllActors,
        e::Exception;
        freeze::Bool=true,
        format::Symbol=:json)

Save FUSE dd, ini, act files and exception stacktrace
"""
function save(
    savedir::AbstractString,
    dd::IMAS.dd,
    ini::ParametersAllInits,
    act::ParametersAllActors,
    e::Exception;
    freeze::Bool=true,
    format::Symbol=:json)

    save(savedir, dd, ini, act; freeze, format)

    open(joinpath(savedir, "error.txt"), "w") do file
        showerror(file, e, catch_backtrace())
    end

    return savedir
end

"""
    load(savedir::AbstractString)

Returns (dd, ini, act) from files read in a folder

`dd` can be in in JSON `dd.json` or HDF `dd.h5` format.

Returns `missing` for files are not there or if `error.txt` file exists in the folder
"""
function load(savedir::AbstractString)
    if isfile(joinpath(savedir, "error.txt"))
        println(savedir)
        return missing, missing, missing
    end
    if isfile(joinpath(savedir, "dd.h5"))
        dd = IMAS.hdf2imas(joinpath(savedir, "dd.h5"))
    elseif isfile(joinpath(savedir, "dd.json"))
        dd = IMAS.json2imas(joinpath(savedir, "dd.json"))
    else
        dd = missing
    end
    if isfile(joinpath(savedir, "ini.json"))
        ini = json2ini(joinpath(savedir, "ini.json"))
    else
        ini = missing
    end
    if isfile(joinpath(savedir, "act.json"))
        act = json2act(joinpath(savedir, "act.json"))
    else
        act = missing
    end
    return dd, ini, act
end

"""
    load(dir::AbstractString, extract::AbstractDict{Symbol,Function})::Dict{Symbol,Any}

Read dd, ini, act from JSON files in a folder and extract some data from them

`extract` functions should accept `dd, ini, act` as inputs, like this:

    extract = Dict(
            :beta_normal => (dd,ini,act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
            :time => (dd,ini,act) -> @ddtime(dd.equilibrium.time),
            :R0 => (dd,ini,act) -> ini.equilibrium.R0
        )

"""
function load(dir::AbstractString, extract::AbstractDict{Symbol,Function})::Dict{Symbol,Any}
    return load(dir, [extract])[1]
end

function load(dir::AbstractString, extracts::AbstractVector{<:AbstractDict{Symbol,Function}})::Vector{Dict{Symbol,Any}}
    dd, ini, act = FUSE.load(dir)
    out = Dict{Symbol,Any}[]
    for extract in extracts
        results = Dict{Symbol,Any}()
        for key in keys(extract)
            if typeof(dd) === typeof(ini) === typeof(act) === Missing
                results[key] = NaN
                continue
            end
            try
                results[key] = extract[key](dd, ini, act)
            catch e
                results[key] = NaN
            end
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
function load(dirs::AbstractVector{<:AbstractString}, extract::Dict{Symbol,Function}; filter_invalid::Bool=true)::DataFrames.DataFrame
    return load(dirs, [extract]; filter_invalid)[1]
end

function load(dirs::AbstractVector{<:AbstractString}, extracts::AbstractVector{<:AbstractDict{Symbol,Function}}; filter_invalid::Bool=true)::Vector{DataFrames.DataFrame}
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
    Threads.@threads for k in eachindex(dirs)
        df[k, :] = load(dirs[k], all_extracts)
        ProgressMeter.next!(p)
    end

    # filter
    if filter_invalid
        df = filter(row -> all(x -> !(x isa Number && (isnan(x) || isinf(x))), row), df)
    end

    # split into separate dataframes
    dfs = DataFrames.DataFrame[]
    for extract in extracts
        push!(dfs, df[:, collect(keys(extract))])
    end

    return dfs
end
