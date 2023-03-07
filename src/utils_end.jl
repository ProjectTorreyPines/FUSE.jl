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

Read dd, ini, act from JSON/HDF files in a folder and extract some data from them

Each of the `extract` functions should accept `dd, ini, act` as inputs, like this:

    extract = Dict(
            :beta_normal => (dd,ini,act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
            :time => (dd,ini,act) -> @ddtime(dd.equilibrium.time),
            :R0 => (dd,ini,act) -> ini.equilibrium.R0
        )

"""
function load(dir::AbstractString, extract::AbstractDict{Symbol,Function})::Dict{Symbol,Any}
    dd, ini, act = FUSE.load(dir)
    results = Dict{Symbol,Any}()
    for key in keys(extract)
        if typeof(dd) === typeof(ini) === typeof(act) === Missing
            results[key] = NaN
            continue
        end
        try
            tmp = extract[key](dd, ini, act)
            if typeof(tmp) <: AbstractDict
                for (k, v) in tmp
                    if typeof(v) <: IMAS.DDigestField
                        v = v.value
                    end
                    results[k] = v
                end
            else
                results[key] = tmp
            end
        catch e
            results[key] = NaN
        end
    end
    return results
end

"""
    load(dirs::AbstractVector{<:AbstractString}, extract::Dict{Symbol,Function})::DataFrames.DataFrame

Read dd, ini, act from JSON/HDF files in multiple directores and extract some data from them returning results in DataFrame format

Each of the `extract` functions should accept `dd, ini, act` as inputs, like this:

    extract = Dict(
            :beta_normal => (dd,ini,act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
            :time => (dd,ini,act) -> @ddtime(dd.equilibrium.time),
            :R0 => (dd,ini,act) -> ini.equilibrium.R0
        )
"""
function load(dirs::AbstractVector{<:AbstractString}, extract::Dict{Symbol,Function}; filter_invalid::Bool=true)::DataFrames.DataFrame
    # allocate memory
    df = DataFrames.DataFrame(load(dirs[1], extract))
    for k in 2:length(dirs)
        push!(df, df[1, :])
    end

    # load the data
    p = ProgressMeter.Progress(length(dirs); showspeed=true)
    Threads.@threads for k in eachindex(dirs)
        df[k, :] = load(dirs[k], extract)
        ProgressMeter.next!(p)
    end

    # filter
    if filter_invalid
        df = filter(row -> all(x -> !(x isa Number && (isnan(x) || isinf(x))), row), df)
    end

    return df
end

"""
    IMAS.digest(dirs::AbstractVector{<:AbstractString}; extract::Dict{Symbol,Function}, filter_invalid::Bool=true)::DataFrames.DataFrame

Digest dd from JSON/HDF files in multiple directores and return results in DataFrame format
"""
function IMAS.digest(dirs::AbstractVector{<:AbstractString}; extract::Dict{Symbol,Function}=Dict{Symbol,Function}(), filter_invalid::Bool=true)::DataFrames.DataFrame
    all_extract = Dict{Symbol,Function}(
        :digest => (dd, ini, act) -> IMAS.digest(dd.summary),
    )
    merge!(all_extract, extract)
    return load(dirs, all_extract; filter_invalid)
end