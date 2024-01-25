using Weave: Weave
using DelimitedFiles: DelimitedFiles

# ========== #
# Checkpoint #
# ========== #
"""
Provides handy checkpoint capabilities of `dd`, `ini`, `act`

It essentially behaves like a ordered dictionary, with the difference that it does deepcopy on getindex and setindex!

Sample usage in a Jupyter notebook:

    chk = FUSE.Checkpoint()
    ...init...
    chk[:init] = dd, ini, act; # store

    --------

    dd, ini, act = chk[:init] # restore
    ...something...
    chk[:something] = dd, ini, act; # store

    --------

    dd = chk[:something].dd # restore only dd
    ...something_else...
    chk[:something_else] = dd # store only dd

    --------

    dd = chk[:something_else].dd # restore only dd
    ...something_else_more...
"""
Base.@kwdef struct Checkpoint
    history::OrderedCollections.OrderedDict = OrderedCollections.OrderedDict()
end

function Base.setindex!(chk::Checkpoint, dd::IMAS.dd, key::Symbol)
    return chk.history[key] = (dd=deepcopy(dd),)
end

function Base.setindex!(chk::Checkpoint, dd_ini_act::Tuple{IMAS.dd{D},ParametersInits{P},ParametersActors{P}}, key::Symbol) where {D<:Real,P<:Real}
    return chk.history[key] = (dd=deepcopy(dd_ini_act[1]), ini=deepcopy(dd_ini_act[2]), act=deepcopy(dd_ini_act[3]))
end

function Base.getindex(chk::Checkpoint, key::Symbol)
    return deepcopy(chk.history[key])
end

function Base.show(io::IO, chk::Checkpoint)
    for (k, v) in chk.history
        TPs = map(typeof, v)
        what = String[]
        if any(tp <: IMAS.dd for tp in TPs)
            push!(what, "dd")
        end
        if any(tp <: ParametersInits for tp in TPs)
            push!(what, "ini")
        end
        if any(tp <: ParametersActors for tp in TPs)
            push!(what, "act")
        end
        println(io, "$(repr(k)) => $(join(what,", "))")
    end
    return
end

# Generate the delegated methods (NOTE: these methods do not require deepcopy)
for func in [:empty!, :delete!, :haskey, :pop!, :popfirst!]
    @eval function Base.$func(chk::Checkpoint, args...; kw...)
        return $func(chk.history, args...; kw...)
    end
end

# ===================================== #
# extract data from FUSE save folder(s) #
# ===================================== #
"""
    IMAS.extract(dir::AbstractString, xtract::T=IMAS.ExtractFunctionsLibrary)::T where {T<:AbstractDict{Symbol,IMAS.ExtractFunction}}

Read dd.json/h5 in a folder and extract data from it.
"""
function IMAS.extract(dir::AbstractString, xtract::T=IMAS.ExtractFunctionsLibrary)::T where {T<:AbstractDict{Symbol,IMAS.ExtractFunction}}
    dd, ini, act = load(dir; load_ini=false, load_act=false, skip_on_error=true)
    return extract(dd, xtract)
end

"""
    IMAS.extract(
        DD::Union{Nothing,Vector{<:Union{AbstractString, IMAS.dd}}},
        xtract::AbstractDict{Symbol, IMAS.ExtractFunction} = IMAS.ExtractFunctionsLibrary;
        filter_invalid::Symbol = :none,
        cache::AbstractString = "",
        read_cache::Bool = true,
        write_cache::Bool = true)::DataFrames.DataFrame

Extract data from multiple FUSE results folders or `dd`s and return results in DataFrame format.

Filtering can by done by `:cols` that have all NaNs, `:rows` that have any NaN, both with `:all`, or `:none`.

Specifying a `cache` file allows caching of extraction results and not having to parse data.

if DD is nothing and cache file is specified, then data is loaded from cachefile alone.
"""
function IMAS.extract(
    DD::Union{Nothing,Vector{<:Union{AbstractString,IMAS.dd}}},
    xtract::AbstractDict{Symbol,IMAS.ExtractFunction}=IMAS.ExtractFunctionsLibrary;
    filter_invalid::Symbol=:none,
    cache::AbstractString="",
    read_cache::Bool=true,
    write_cache::Bool=true)::DataFrames.DataFrame

    # test filter_invalid
    @assert filter_invalid ∈ (:none, :cols, :rows, :all) "filter_invalid can only be one of [:none, :cols, :rows, :all]"

    if DD !== nothing && length(cache) > 0
        @assert typeof(DD[1]) <: AbstractString "cache is only meant to work when extracting data from FUSE results folders"
    end

    # load in cache
    cached_dirs = []
    if length(cache) > 0 && read_cache && isfile(cache)
        df_cache = DataFrames.DataFrame(CSV.File(cache))
        cached_dirs = df_cache[:, :dir]
        @info "Loaded cache file with $(length(cached_dirs)) results"
    end

    if DD === nothing
        df = df_cache

    else
        # allocate memory
        tmp = Dict(extract(DD[1], xtract))
        tmp[:dir] = abspath(DD[1])
        df = DataFrames.DataFrame(tmp)
        for k in 2:length(DD)
            push!(df, df[1, :])
        end

        # load the data
        p = ProgressMeter.Progress(length(DD); showspeed=true)
        Threads.@threads for k in eachindex(DD)
            aDDk = abspath(DD[k])
            try
                if aDDk in cached_dirs
                    kcashe = findfirst(dir -> dir == aDDk, cached_dirs)
                    df[k, :] = df_cache[kcashe, :]
                else
                    tmp = Dict(extract(aDDk, xtract))
                    tmp[:dir] = aDDk
                    df[k, :] = tmp
                end
            catch
                continue
            end
            ProgressMeter.next!(p)
        end
        ProgressMeter.finish!(p)

        # cache extrated information to CSV file
        if length(cache) > 0 && write_cache
            DelimitedFiles.writedlm(cache, Iterators.flatten(([names(df)], eachrow(df))), ',')
            @info "Written cache file with $(length(df.dir)) results"
        end
    end

    # filter
    if filter_invalid ∈ (:cols, :all)
        # drop columns that have all NaNs
        isnan_nostring(x::Any) = (typeof(x) <: Number) ? isnan(x) : false
        visnan(x::AbstractVector) = isnan_nostring.(x)
        df = df[:, .!all.(visnan.(eachcol(df)))]
    end
    if filter_invalid ∈ (:rows, :all)
        # drop rows that have any NaNs
        df = filter(row -> all(x -> !(x isa Number && (isnan(x) || isinf(x))), row), df)
    end

    return df
end

"""
    DataFrames.DataFrame(xtract::AbstractDict{Symbol,IMAS.ExtractFunction})

Construct a DataFrame from a dictionary of IMAS.ExtractFunction
"""
function DataFrames.DataFrame(xtract::AbstractDict{Symbol,IMAS.ExtractFunction})
    return DataFrames.DataFrame(Dict(xtract))
end

"""
    Dict(xtract::AbstractDict{Symbol,IMAS.ExtractFunction})

Construct a Dictionary with the evaluated values of a dictionary of IMAS.ExtractFunction
"""
function Dict(xtract::AbstractDict{Symbol,IMAS.ExtractFunction})
    tmp = Dict()
    for xfun in values(xtract)
        tmp[xfun.name] = xfun.value
    end
    return tmp
end

# ==================== #
# save/load simulation #
# ==================== #
"""
    save(
        savedir::AbstractString,
        dd::Union{Nothing,IMAS.dd},
        ini::Union{Nothing,ParametersAllInits},
        act::Union{Nothing,ParametersAllActors},
        e::Union{Nothing,Exception}=nothing;
        timer::Bool=true,
        varinfo::Bool=false,
        freeze::Bool=true,
        format::Symbol=:json,
        overwrite_files::Bool=true)

Save FUSE (`dd`, `ini`, `act`) to `dd.json`/`h5`, `ini.json`, and `act.json` files and exception stacktrace to `error.txt`

`timer` option allows saving of `FUSE.timer` info to `timer.txt` file

`varinfo` option allows saving of detailed variables memory usage to `varinfo.txt` file

If `dd`, `ini`, `act`, or `e` are `nothing` then the corresponding file is not created.
"""
function save(
    savedir::AbstractString,
    dd::Union{Nothing,IMAS.dd},
    ini::Union{Nothing,ParametersAllInits},
    act::Union{Nothing,ParametersAllActors},
    e::Union{Nothing,Exception}=nothing;
    timer::Bool=true,
    varinfo::Bool=false,
    memtrace::Bool=false,
    freeze::Bool=true,
    format::Symbol=:json,
    overwrite_files::Bool=true)

    @assert format ∈ (:hdf, :json) "format must be either `:hdf` or `:json`"

    savedir = abspath(savedir)
    if !overwrite_files || !isdir(savedir)
        mkdir(savedir)
    end

    # first write error.txt so that if we are parsing while running optimizer,
    # the parser can immediately see if this is a failing case
    if e !== nothing
        open(joinpath(savedir, "error.txt"), "w") do file
            return showerror(file, e, catch_backtrace())
        end
    end

    # save timer output
    if timer
        open(joinpath(savedir, "timer.txt"), "w") do file
            return show(file, FUSE.timer)
        end
    end

    # save vars usage
    if varinfo
        open(joinpath(savedir, "varinfo.txt"), "w") do file
            return println(file, FUSE.varinfo(FUSE; all=true, imported=true, recursive=true, sortby=:size, minsize=1024))
        end
    end

    # save memory trace
    if memtrace
        save(FUSE.memtrace, joinpath(savedir, "memtrace.txt"))
    end

    if dd !== nothing
        if format == :hdf
            IMAS.imas2hdf(dd, joinpath(savedir, "dd.h5"); freeze)
        elseif format == :json
            IMAS.imas2json(dd, joinpath(savedir, "dd.json"); freeze)
        end
    end

    if ini !== nothing
        ini2json(ini, joinpath(savedir, "ini.json"))
    end

    if act !== nothing
        act2json(act, joinpath(savedir, "act.json"))
    end

    return savedir
end

"""
    load(savedir::AbstractString; load_dd::Bool=true, load_ini::Bool=true, load_act::Bool=true, skip_on_error::Bool=false)

Read (dd, ini, act) to dd.json/h5, ini.json, and act.json files.

Returns `missing` for files are not there or if `error.txt` file exists in the folder.
"""
function load(savedir::AbstractString; load_dd::Bool=true, load_ini::Bool=true, load_act::Bool=true, skip_on_error::Bool=false)
    if skip_on_error && isfile(joinpath(savedir, "error.txt"))
        @warn "$savedir simulation errored"
        return missing, missing, missing
    end
    dd = missing
    if load_dd
        if isfile(joinpath(savedir, "dd.h5"))
            dd = IMAS.hdf2imas(joinpath(savedir, "dd.h5"))
        elseif isfile(joinpath(savedir, "dd.json"))
            dd = IMAS.json2imas(joinpath(savedir, "dd.json"))
        end
    end
    ini = missing
    if load_ini && isfile(joinpath(savedir, "ini.json"))
        ini = json2ini(joinpath(savedir, "ini.json"))
    end
    act = missing
    if load_act && isfile(joinpath(savedir, "act.json"))
        act = json2act(joinpath(savedir, "act.json"))
    end
    return dd, ini, act
end

"""
    digest(
        dd::IMAS.dd;
        terminal_width::Int=136,
        line_char::Char='─',
        section::Int=0)

Provides concise and informative summary of `dd`, including several plots

NOTES:

  - `section` is used internally to produce digest PDFs
  - this function is defined in FUSE and not IMAS because it uses Plots.jl and not BaseRecipies.jl
"""
function digest(
    dd::IMAS.dd;
    terminal_width::Int=136,
    line_char::Char='─',
    section::Int=0)

    sec = 1
    if section ∈ (0, sec)
        IMAS.print_tiled(extract(dd); terminal_width, line_char)
    end

    # equilibrium with build and PFs
    sec += 1
    if !isempty(dd.equilibrium.time_slice) && section ∈ (0, sec)
        println('\u200B')
        p = plot(dd.equilibrium; legend=false)
        if !isempty(dd.build.layer)
            plot!(p[1], dd.build; legend=false)
        end
        if !isempty(dd.pf_active.coil)
            plot!(p[1], dd.pf_active; legend=false, colorbar=false)
        end
        if !isempty(dd.divertors.divertor)
            plot!(p[1], dd.divertors; legend=false)
        end
        display(p)
    end

    # build layers
    sec += 1
    if !isempty(dd.build.layer) && section ∈ (0, sec)
        println('\u200B')
        display(dd.build.layer)
    end

    # core profiles
    # temperatures
    sec += 1
    if !isempty(dd.core_profiles.profiles_1d) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_profiles; only=1))
    end
    # densities
    sec += 1
    if !isempty(dd.core_profiles.profiles_1d) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_profiles; only=2))
    end
    # rotation
    sec += 1
    if !isempty(dd.core_profiles.profiles_1d) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_profiles; only=3))
    end

    # core sources
    # electron heat
    sec += 1
    if !isempty(dd.core_sources.source) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_sources; only=1))
    end
    # ion heat
    sec += 1
    if !isempty(dd.core_sources.source) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_sources; only=2))
    end
    # electron particle
    sec += 1
    if !isempty(dd.core_sources.source) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_sources; only=3))
    end
    # parallel current
    sec += 1
    if !isempty(dd.core_sources.source) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_sources; only=4))
    end

    # Electron energy flux matching
    sec += 1
    if !isempty(dd.core_transport) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_transport; only=1))
    end

    # Ion energy flux matching
    sec += 1
    if !isempty(dd.core_transport) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_transport; only=2))
    end

    # Electron particle flux matching
    sec += 1
    if !isempty(dd.core_transport) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_transport; only=3))
    end

    # Momentum flux matching
    sec += 1
    if !isempty(dd.core_transport) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.core_transport; only=4))
    end

    # neutron wall loading
    sec += 1
    if !isempty(dd.neutronics.time_slice) && section ∈ (0, sec)
        println('\u200B')
        xlim = extrema(dd.neutronics.first_wall.r)
        xlim = (xlim[1] - ((xlim[2] - xlim[1]) / 10.0), xlim[2] + ((xlim[2] - xlim[1]) / 10.0))
        l = @layout [a{0.3w} b{0.6w,0.9h}]
        p = plot(; layout=l, size=(900, 400))
        plot!(p, dd.neutronics.time_slice[].wall_loading; xlim, subplot=1)
        neutrons = define_neutrons(dd, 100000)[1]
        plot!(p, neutrons, dd.equilibrium.time_slice[]; xlim, subplot=1)
        plot!(p, dd.neutronics.time_slice[].wall_loading; cx=false, subplot=2, ylabel="")
        display(p)
    end

    # SOL
    sec += 1
    if !isempty(dd.equilibrium.time_slice) && section ∈ (0, sec)
        println('\u200B')
        p = plot(dd.wall)
        plot!(IMAS.sol(dd); xlim=[0.0, Inf])
        display(p)
    end

    # center stack stresses
    sec += 1
    if !ismissing(dd.solid_mechanics.center_stack.grid, :r_oh) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.solid_mechanics.center_stack.stress))
    end

    # pf active
    sec += 1
    if !isempty(dd.pf_active.coil) && section ∈ (0, sec)
        println('\u200B')
        time0 = dd.equilibrium.time[end]
        l = @layout [a{0.5w} b{0.5w}]
        p = plot(; layout=l, size=(900, 400))
        plot!(p, dd.pf_active, :currents; time0, title="PF currents at t=$(time0) s", subplot=1)
        plot!(p, dd.equilibrium; time0, cx=true, subplot=2)
        plot!(p, dd.build; subplot=2, legend=false)
        plot!(p, dd.pf_active; time0, subplot=2)
        display(p)
    end

    # tf
    sec += 1
    if !ismissing(dd.build.tf, :coils_n) && !isempty(dd.build.layer) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.build.tf; title="TF coils -- Top view"))
    end

    # balance of plant (cannot be plotted right now plotting can only be done when running actor and not from data in dd)
    # sec += 1
    # if !missing(dd.balance_of_plant, :Q_plant) && section ∈ (0, sec)
    # println('\u200B')
    #     display(plot(dd.balance_of_plant))
    # end

    # costing
    sec += 1
    if !ismissing(dd.costing.cost_direct_capital, :cost) && (dd.costing.cost_direct_capital.cost != 0) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.costing.cost_direct_capital))
    end

    return nothing
end

"""
    digest(dd::IMAS.dd,
        title::AbstractString,
        description::AbstractString="";
        ini::Union{Nothing,ParametersAllInits}=nothing,
        act::Union{Nothing,ParametersAllActors}=nothing)

Write digest to PDF in current working directory.

PDF filename is based on title (with `" "` replaced by `"_"`)
"""
function digest(dd::IMAS.dd,
    title::AbstractString,
    description::AbstractString="";
    ini::Union{Nothing,ParametersAllInits}=nothing,
    act::Union{Nothing,ParametersAllActors}=nothing
)
    title = replace(title, r".pdf$" => "", "_" => " ")
    outfilename = joinpath(pwd(), "$(replace(title," "=>"_")).pdf")

    tmpdir = mktempdir()
    logger = SimpleLogger(stderr, Logging.Warn)
    try
        filename = redirect_stdout(Base.DevNull()) do
            filename = with_logger(logger) do
                return Weave.weave(joinpath(@__DIR__, "digest.jmd");
                    latex_cmd=["xelatex"],
                    mod=@__MODULE__,
                    doctype="md2pdf",
                    template=joinpath(@__DIR__, "digest.tpl"),
                    out_path=tmpdir,
                    args=Dict(
                        :dd => dd,
                        :ini => ini,
                        :act => act,
                        :title => title,
                        :description => description))
            end
        end
        cp(filename, outfilename; force=true)
        return outfilename
    catch e
        println("Generation of $(basename(outfilename)) failed. See directory: $tmpdir\n$e")
    else
        rm(tmpdir; recursive=true, force=true)
    end
end

"""
    categorize_errors(dirs::AbstractVector{<:AbstractString}; show_first_line=false, do_plot=true, extra_error_messages::AbstractDict=Dict())

Looks at the first line of each error.txt file in dirs and categorizes them
"""
function categorize_errors(
    dirs::AbstractVector{<:AbstractString};
    show_first_line=false,
    do_plot=true,
    extra_error_messages::AbstractDict{<:AbstractString,Symbol}=Dict{String,Symbol}()
)
    # error counting and error message dict
    errors = Dict(:other => String[])
    error_messages = Dict(
        "EQDSK_COCOS_01.OUT" => :chease,
        "Unable to blend the core-pedestal" => :blend_core_ped,
        "Bad expression" => :bad_expression,
        "Exceeded limits" => :exceed_lim_A,
        "Some stability models have breached their limit threshold:" => :exceed_lim_B,
        "TaskFailedException" => :task_exception,
        "Could not trace closed flux surface" => :flux_surfaces_A,
        "Flux surface at ψ=" => :flux_surfaces_B,
        "OH stresses" => :OH_stresses,
        "TF stresses" => :TF_stresses,
        "yield_strength" => :CS_stresses,
        "TF cannot achieve requested B0" => :TF_limit,
        "The OH flux is insufficient to have any flattop duration" => :OH_flux,
        "OH cannot achieve requested flattop" => :OH_flattop,
        "OH exceeds critical current" => :OH_critical_j,
        "< dd.build.tf.critical_j" => :TF_critical_j,
        "DomainError with" => :Solovev,
        "BoundsError: attempt to access" => :flux_surfaces_C,
        "divertors" => :divertors)
    merge!(error_messages, extra_error_messages)

    other_errors = Dict{String,Vector{String}}()

    # go through directories
    for dir in dirs
        filename = joinpath([dir, "error.txt"])
        if !isfile(filename)
            continue
        end
        first_line, second_line = open(filename, "r") do f
            return (readline(f), readline(f))
        end
        found = false
        for (err, cat) in error_messages
            if occursin(err, first_line)
                found = true
                if cat == :exceed_lim_B
                    cat = Symbol(second_line)
                end
                if cat ∉ keys(errors)
                    errors[cat] = String[]
                end
                push!(errors[cat], dir)
                break
            end
        end
        if !found
            push!(errors[:other], dir)
            if "$(first_line)\n$(second_line)" ∉ keys(other_errors)
                other_errors["$(first_line)\n$(second_line)"] = String[]
            end
            push!(other_errors["$(first_line)\n$(second_line)"], dir)
        end
    end

    if do_plot
        display(histogram(findall(x -> isfile(joinpath(x, "error.txt")), sort(dirs)); labe="Errors"))
        labels = collect(keys(errors))
        v = collect(map(length, values(errors)))
        index = sortperm(v)[end:-1:1]
        display(pie(["$(rpad(string(length(errors[cat])),8))   $(string(cat))" for cat in labels[index]], v[index]; legend=:outerright))
    end

    if show_first_line
        for (k, v) in other_errors
            println(k)
            println(v)
            println()
        end
    end

    return errors
end
