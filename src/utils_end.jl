
mutable struct ExtractFunction
    group::Symbol
    name::Symbol
    units::String
    func::Function
    # inner constructor to register ExtractFunction in ExtractFunctionsLibrary
    ExtractFunction(group::Symbol, name::Symbol, units::String, func::Function) = begin
        objf = new(group, name, units, func)
        ExtractFunctionsLibrary[objf.name] = objf
        return objf
    end
end

const ExtractFunctionsLibrary = DataStructures.OrderedDict{Symbol,ExtractFunction}()
function update_ExtractFunctionsLibrary!()
    empty!(ExtractFunctionsLibrary)
    ExtractFunction(:equilibrium, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
    ExtractFunction(:equilibrium, :δ, "-", dd -> dd.equilibrium.time_slice[].boundary.triangularity)
    ExtractFunction(:equilibrium, :ζ, "-", dd -> dd.equilibrium.time_slice[].boundary.squareness)
    ExtractFunction(:equilibrium, :B0, "T", dd -> @ddtime(dd.summary.global_quantities.b0.value))
    ExtractFunction(:equilibrium, :ip, "MA", dd -> @ddtime(dd.summary.global_quantities.ip.value)/1e6)
    ExtractFunction(:equilibrium, :R0, "m", dd -> dd.summary.global_quantities.r0.value)
    ExtractFunction(:equilibrium, :βn, "-", dd -> @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
    ExtractFunction(:profiles, :zeff, "-", dd -> @ddtime(dd.summary.volume_average.zeff.value))
    ExtractFunction(:profiles, :Te0, "keV", dd -> dd.core_profiles.profiles_1d[].electrons.temperature[1] / 1E3)
    ExtractFunction(:profiles, :Ti0, "keV", dd -> dd.core_profiles.profiles_1d[].ion[1].temperature[1] / 1E3)
    ExtractFunction(:profiles, :Pfusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6)
    ExtractFunction(:balance_of_plant, :Pelectric_net, "MWe", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6)
    ExtractFunction(:heating_current_drive, :Pelectron_cyclotron, "W", dd -> @ddtime(summary.heating_current_drive.power_launched_ec))
    ExtractFunction(:heating_current_drive, :Pneutral_beam, "W", dd -> @ddtime(summary.heating_current_drive.power_launched_nbi))
    ExtractFunction(:heating_current_drive, :Pion_cyclotron, "W", dd -> @ddtime(summary.heating_current_drive.power_launched_ic))
    ExtractFunction(:heating_current_drive, :Plower_hybrid, "W", dd -> @ddtime(summary.heating_current_drive.power_launched_lh))
    ExtractFunction(:heating_current_drive, :Paux_total, "W", dd -> @ddtime(summary.heating_current_drive.power_launched_total))
    ExtractFunction(:costing, :levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE)
    ExtractFunction(:costing, :capital_cost, "\$M", dd -> dd.costing.cost_direct_capital.cost)
    ExtractFunction(:build, :flattop, "\$M", dd -> dd.build.oh.flattop_duration / 3600.0)
end
update_ExtractFunctionsLibrary!()

"""
    (ef::ExtractFunction)(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

run the extract function
"""
function (ef::ExtractFunction)(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    return ef.func(dd)
end

"""
    (ef::ExtractFunction)(dd::IMAS.dd)

run the extract function
"""
function (ef::ExtractFunction)(dd::IMAS.dd)
    return ef.func(dd)
end

function Base.show(io::IO, f::ExtractFunction)
    printstyled(io, f.group; bold=true)
    printstyled(io, "."; bold=true)
    printstyled(io, f.name; bold=true, color=:blue)
    print(io, " →")
end

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
    load(dir::AbstractString, extract::AbstractDict{Symbol,T})::Dict{Symbol,Any} where {T<:Union{Function,ExtractFunction}}

Read dd, ini, act from JSON/HDF files in a folder and extract some data from them

Each of the `extract` functions should accept `dd, ini, act` as inputs, like this:

    extract = Dict(
            :beta_normal => (dd,ini,act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
            :time => (dd,ini,act) -> @ddtime(dd.equilibrium.time),
            :R0 => (dd,ini,act) -> ini.equilibrium.R0
        )

"""
function load(dir::AbstractString, extract::AbstractDict{Symbol,T})::Dict{Symbol,Any} where {T<:Union{Function,ExtractFunction}}
    dd, ini, act = FUSE.load(dir)
    results = Dict{Symbol,Any}()
    for key in keys(extract)
        if typeof(dd) === typeof(ini) === typeof(act) === Missing
            results[key] = NaN
            continue
        end
        # try
            results[key] = extract[key](dd, ini, act)
        # catch e
        #     rethrow(e)
        #     results[key] = NaN
        # end
    end
    return results
end

"""
    load(dirs::AbstractVector{<:AbstractString}, extract::AbstractDict{Symbol,T}; filter_invalid::Bool=true)::DataFrames.DataFrame where {T<:Union{Function,ExtractFunction}}

Read dd, ini, act from JSON/HDF files in multiple directores and extract some data from them returning results in DataFrame format

Each of the `extract` functions should accept `dd, ini, act` as inputs, like this:

    extract = Dict(
            :beta_normal => (dd,ini,act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
            :time => (dd,ini,act) -> @ddtime(dd.equilibrium.time),
            :R0 => (dd,ini,act) -> ini.equilibrium.R0
        )
"""
function load(dirs::AbstractVector{<:AbstractString}, extract::AbstractDict{Symbol,T}=ExtractFunctionsLibrary; filter_invalid::Bool=true)::DataFrames.DataFrame where {T<:Union{Function,ExtractFunction}}
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
