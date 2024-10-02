import Weave
using InteractiveUtils: summarysize, format_bytes, Markdown
import DelimitedFiles
import OrderedCollections
import DataFrames
import ProgressMeter

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

also works with @checkin and @checkout macros

    chk = FUSE.Checkpoint()
    ...init...
    @checkin chk :init dd ini act # store dd, ini, act

    --------

    @checkout chk :init dd ini act # restore dd, ini, act
    ...something...
    @checkin chk :something dd act # store dd, act

    --------

    @checkout chk :something dd # restore only dd
    ...something_else...
    @checkin chk :something_else dd # store only dd

    --------

    @checkout chk :something_else dd # restore only dd
    ...something_else_more...
"""
Base.@kwdef struct Checkpoint
    history::OrderedCollections.OrderedDict = OrderedCollections.OrderedDict()
end

function Base.setindex!(chk::Checkpoint, dd::IMAS.dd, key::Symbol)
    return chk.history[key] = (dd=deepcopy(dd),)
end

function Base.setindex!(chk::Checkpoint, dd_ini_act::Tuple{IMAS.dd,ParametersAllActors}, key::Symbol)
    return chk.history[key] = (dd=deepcopy(dd_ini_act[1]), act=deepcopy(dd_ini_act[2]))
end

function Base.setindex!(chk::Checkpoint, dd_ini_act::Tuple{IMAS.dd,ParametersAllInits,ParametersAllActors}, key::Symbol)
    return chk.history[key] = (dd=deepcopy(dd_ini_act[1]), ini=deepcopy(dd_ini_act[2]), act=deepcopy(dd_ini_act[3]))
end

function Base.getindex(chk::Checkpoint, key::Symbol)
    return deepcopy(chk.history[key])
end

function Base.show(io::IO, ::MIME"text/plain", chk::Checkpoint)
    for (k, v) in chk.history
        what = Tuple{Symbol,Type}[(k, typeof(v)) for (k, v) in pairs(v)]
        println(io, "$(repr(k)) => $(join(what,", "))")
    end
    return nothing
end

# Generate the delegated methods (NOTE: these methods do not require deepcopy)
for func in [:empty!, :delete!, :haskey, :pop!, :popfirst!]
    @eval function Base.$func(chk::Checkpoint, args...; kw...)
        return $func(chk.history, args...; kw...)
    end
end

"""
    @checkin chk :key a b c

Macro to save variables into a Checkpoint under a specific key
"""
macro checkin(checkpoint, key, vars...)
    key = esc(key)
    checkpoint = esc(checkpoint)

    # Save all the variables in the `vars` list under the provided key using their names
    return quote
        d = getfield($checkpoint, :history)
        if $key in keys(d)
            dict = Dict(k => v for (k, v) in pairs(d[$key]))
        else
            dict = Dict{Symbol,Any}()
        end
        $(Expr(:block, [:(dict[Symbol($(string(v)))] = deepcopy($(esc(v)))) for v in vars]...))
        d[$key] = NamedTuple{Tuple(keys(dict))}(values(dict))  # Convert the dictionary to a NamedTuple
        nothing
    end
end

"""
    @checkout chk :key a c

Macro to load variables from a Checkpoint
"""
macro checkout(checkpoint, key, vars...)
    key = esc(key)
    checkpoint = esc(checkpoint)

    # Restore variables from the checkpoint
    return quote
        d = getfield($checkpoint, :history)
        if haskey(d, $key)
            saved_vars = d[$key]
            $(Expr(:block, [:($(esc(v)) = deepcopy(getfield(saved_vars, Symbol($(string(v)))))) for v in vars]...))
        else
            throw(KeyError($key))
        end
        nothing
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
    IMAS.last_time(dd)
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
Filtering by `:cols_row1` removes the columns that have NaN in the first row of the table.

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

    function identifier(DD::AbstractVector{<:AbstractString}, k::Int)
        return abspath(DD[k]), abspath(DD[k])
    end

    function identifier(DD::AbstractVector{<:IMAS.dd}, k::Int)
        return "$k", DD[k]
    end

    # test filter_invalid
    @assert filter_invalid ∈ (:none, :cols_row1, :cols, :rows, :all) "filter_invalid can only be one of [:none, :cols_row1, :cols, :rows, :all]"

    if DD !== nothing && length(cache) > 0
        @assert typeof(DD[1]) <: AbstractString "cache is only meant to work when extracting data from FUSE results folders"
    end

    # load in cache
    cached_dirs = []
    df_cache = DataFrames.DataFrame()
    if length(cache) > 0 && read_cache && isfile(cache)
        df_cache = DataFrames.DataFrame(CSV.File(cache))
        cached_dirs = df_cache[:, :dir]
        @info "Loaded cache file with $(length(cached_dirs)) results"
    end

    if DD === nothing
        df = df_cache

    else
        # allocate memory
        if filter_invalid == :cols_row1
            xtract = deepcopy(xtract)
            xtr = extract(DD[1], xtract)
            tmp = Dict()
            for (xkey, xfun) in xtr
                if xfun.value === NaN
                    delete!(xtract, xkey)
                else
                    tmp[xfun.name] = xfun.value
                end
            end
        else
            tmp = Dict(extract(DD[1], xtract))
        end

        tmp[:dir] = ""
        df = DataFrames.DataFrame(tmp)
        for k in 2:length(DD)
            push!(df, df[1, :])
        end

        # load the data
        p = ProgressMeter.Progress(length(DD); showspeed=true)
        Threads.@threads for k in eachindex(DD)
            aDDk, aDD = identifier(DD, k)
            try
                if aDDk in cached_dirs
                    k_cache = findfirst(dir -> dir == aDDk, cached_dirs)
                    df[k, :] = df_cache[k_cache, :]
                else
                    tmp1 = Dict(extract(aDD, xtract))
                    for key in keys(tmp1)
                        if typeof(tmp1[key]) <: Float64 && typeof(tmp[key]) <: String
                            tmp1[key] = ""
                        end
                    end
                    tmp1[:dir] = aDDk
                    df[k, :] = tmp1
                end
            catch e
                if isa(e, InterruptException)
                    rethrow(e)
                end
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
        act::Union{Nothing,ParametersAllActors};
        error::Any=nothing,
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
    act::Union{Nothing,ParametersAllActors};
    error::Any=nothing,
    timer::Bool=true,
    varinfo::Bool=false,
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
    if typeof(error) <: Nothing
        # pass
    elseif typeof(error) <: Exception
        open(joinpath(savedir, "error.txt"), "w") do file
            return showerror(file, error, catch_backtrace())
        end
    else
        open(joinpath(savedir, "error.txt"), "w") do file
            return println(file, string(error))
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
    if parse(Bool, get(ENV, "FUSE_MEMTRACE", "false"))
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

Read (`dd`, `ini`, `act`) from `dd.json/h5`, `ini.json/yaml`, and `act.json/yaml` files.

Returns `missing` for files are not there or if `error.txt` file exists in the folder.
"""
function load(savedir::AbstractString; load_dd::Bool=true, load_ini::Bool=true, load_act::Bool=true, skip_on_error::Bool=false)
    if skip_on_error && isfile(joinpath(savedir, "error.txt"))
        @warn "$savedir simulation errored"
        return missing, missing, missing
    end

    # dd
    if load_dd && isfile(joinpath(savedir, "dd.h5"))
        dd = IMAS.hdf2imas(joinpath(savedir, "dd.h5"))
    elseif load_dd && isfile(joinpath(savedir, "dd.json"))
        dd = IMAS.json2imas(joinpath(savedir, "dd.json"))
    else
        dd = missing
    end

    # ini
    if load_ini && isfile(joinpath(savedir, "ini.json"))
        ini = json2ini(joinpath(savedir, "ini.json"))
    elseif load_ini && isfile(joinpath(savedir, "ini.yaml"))
        ini = yaml2ini(joinpath(savedir, "ini.yaml"))
    else
        ini = missing
    end

    # act
    if load_act && isfile(joinpath(savedir, "act.json"))
        act = json2act(joinpath(savedir, "act.json"))
    elseif load_act && isfile(joinpath(savedir, "act.yaml"))
        act = yaml2act(joinpath(savedir, "act.yaml"))
    else
        act = missing
    end

    return (dd=dd, ini=ini, act=act)
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
        println("@ time = $(dd.global_time) [s]")
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
    for k in 1:5+length(IMAS.list_ions(dd.core_sources))
        if !isempty(dd.core_sources.source) && section ∈ (0, sec)
            println('\u200B')
            display(plot(dd.core_sources; only=k))
        end
    end

    # core sources
    for k in 1:4+length(IMAS.list_ions(dd.core_transport))
        if !isempty(dd.core_transport) && section ∈ (0, sec)
            println('\u200B')
            display(plot(dd.core_transport; only=k))
        end
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
        plot!(p, neutrons, dd.equilibrium.time_slice[]; xlim, subplot=1, colorbar_entry=false)
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
    if !isempty(dd.pf_active.coil) && !ismissing(dd.equilibrium, :time) && section ∈ (0, sec)
        println('\u200B')
        time0 = dd.equilibrium.time[end]
        l = @layout [a{0.5w} b{0.5w}]
        p = plot(; layout=l, size=(900, 400))
        plot!(p, dd.pf_active, :currents; time0, title="PF currents at t=$(time0) s", subplot=1)
        plot!(p, dd.equilibrium; time0, cx=true, subplot=2)
        plot!(p, dd.build; subplot=2, legend=false, equilibrium=false, pf_active=false)
        plot!(p, dd.pf_active; time0, subplot=2, coil_names=true)
        plot!(p, dd.build.pf_active.rail; subplot=2)
        display(p)
    end

    # circuits
    sec += 1
    if !isempty(dd.pf_active.circuit) && section ∈ (0, sec)
        println('\u200B')
        l = @layout length(dd.pf_active.circuit)
        p = plot(; layout=l, size=(900, 900))
        for (k, circuit) in enumerate(dd.pf_active.circuit)
            plot!(p, circuit; subplot=k)
        end
        display(p)
    end

    # pulse_schedule
    sec += 1
    if !isempty(dd.pulse_schedule) && section ∈ (0, sec)
        println('\u200B')
        display(plot(dd.pulse_schedule; title="Pulse schedule"))
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
        if isa(e, InterruptException)
            rethrow(e)
        end
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
        "DomainError with" => :some_negative_root,
        "AssertionError: The output flux is NaN check your transport model fluxes" => :issue_with_transport,
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

# ====== #
# Memory #
# ====== #
# NOTE: Memory tracking is disabled by default, since it has a performance hit
#       To enable, do: export FUSE_MEMTRACE=true"
#
Base.@kwdef struct MemTrace
    data::Vector{Tuple{Dates.DateTime,String,Int}} = Tuple{Dates.DateTime,String,Int}[]
end

const memtrace = MemTrace()

function memory_time_tag(txt::String)
    return push!(memtrace.data, (Dates.now(), txt, get_julia_process_memory_usage()))
end

function memory_time_tag(actor::AbstractActor, msg::String)
    if parse(Bool, get(ENV, "FUSE_MEMTRACE", "false"))
        txt = "$(name(actor)) - @$msg"
        return push!(memtrace.data, (Dates.now(), txt, get_julia_process_memory_usage()))
    end
end

"""
    plot_memtrace(memtrace::MemTrace, n_big_jumps::Int=5, ignore_first_seconds::Int=0)

Plot the memory usage over time from a `MemTrace` object.

# Arguments

  - `n_big_jumps`: number of significant memory jumps to highlight in the plot.
  - `ignore_first_seconds`: number of initial seconds to ignore in the plot. Memory usage will be plotted relative to the memory after this cutoff. Default is `0` (no seconds are ignored).

# Returns

A plot with the following characteristics:

  - Time is displayed on the x-axis as delta seconds since the first recorded time (after ignoring the specified initial seconds).
  - Memory usage is displayed on the y-axis in MB.
  - Memory usage is shown as a scatter plot.
  - The `n_big_jumps` largest jumps in memory are highlighted in red with annotations indicating the action causing each jump.
"""
@recipe function plot_memtrace(memtrace::MemTrace, n_big_jumps::Int=5, ignore_first_seconds::Int=0)
    if isempty(memtrace.data)
        cutoff = 0.0
    else
        cutoff = memtrace.data[1][1] + Dates.Second(ignore_first_seconds)
    end
    filtered_data = filter(point -> point[1] > cutoff, memtrace.data)

    dates = [point[1] for point in filtered_data]
    if !isempty(memtrace.data)
        dates = Dates.value.(dates .- dates[1]) ./ 1000
    end
    action = [point[2] for point in filtered_data]
    mem = [point[3] / 1024 / 1024 for point in filtered_data]
    if !isempty(memtrace.data) && ignore_first_seconds > 0
        mem = mem .- mem[1]
    end

    @series begin
        seriestype := scatter
        label --> ""
        dates, mem
    end

    index = sortperm(diff(mem))[end-min(n_big_jumps - 1, length(mem) - 2):end]
    for i in index
        @series begin
            primary := false
            series_annotations := Plots.text.([action[i], action[i+1]], 6, :red)
            color := :red
            [dates[i], dates[i+1]], [mem[i], mem[i+1]]
        end
    end
    @series begin
        primary := false
        xlabel := "ΔTime [s]"
        ylabel := (ignore_first_seconds > 0 ? "Δ" : "") * "Memory [MB]"
        [], []
    end
end

"""
    get_julia_process_memory_usage()

Returns memory used by current julia process
"""
function get_julia_process_memory_usage()
    pid = getpid()
    mem_info = read(`ps -p $pid -o rss=`, String)
    mem_usage_kb = parse(Int, strip(mem_info))
    return mem_usage_kb * 1024
end

"""
    save(memtrace::MemTrace, filename::String="memtrace.txt")

Save a FUSE memory trace to file
"""
function save(memtrace::MemTrace, filename::String="memtrace.txt")
    open(filename, "w") do file
        for (date, txt, kb) in memtrace.data
            println(file, "$date $kb \"$txt\"")
        end
    end
end

"""
    load(memtrace::MemTrace, filename::String="memtrace.txt")

Load a FUSE memory trace from file
"""
function load(memtrace::MemTrace, filename::String="memtrace.txt")
    open(filename, "r") do file
        for line in eachline(file)
            m = match(r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3})\s+(\d+)\s+\"(.*?)\"", line)
            if m !== nothing
                date = Dates.DateTime(m[1], "yyyy-mm-ddTHH:MM:SS.sss")
                kb = parse(Int, m[2])
                txt = m[3]
                push!(memtrace.data, (date, txt, kb))
            end
        end
    end
    return memtrace
end

"""
    varinfo(m::Module=Main, pattern::Regex=r""; all=false, imported=false, recursive=false, sortby::Symbol=:name, minsize::Int=0)

Return a markdown table giving information about exported global variables in a module, optionally restricted
to those matching `pattern`.

The memory consumption estimate is an approximate lower bound on the size of the internal structure of the object.

  - `all` : also list non-exported objects defined in the module, deprecated objects, and compiler-generated objects.
  - `imported` : also list objects explicitly imported from other modules.
  - `recursive` : recursively include objects in sub-modules, observing the same settings in each.
  - `sortby` : the column to sort results by. Options are `:name` (default), `:size`, and `:summary`.
  - `minsize` : only includes objects with size at least `minsize` bytes. Defaults to `0`.

The output of `varinfo` is intended for display purposes only.  See also [`names`](@ref) to get an array of symbols defined in
a module, which is suitable for more general manipulations.
"""
function varinfo(m::Module=Main, pattern::Regex=r""; all::Bool=false, imported::Bool=false, recursive::Bool=false, sortby::Symbol=:name, minsize::Int=0)
    sortby in (:name, :size, :summary) || throw(ArgumentError("Unrecognized `sortby` value `:$sortby`. Possible options are `:name`, `:size`, and `:summary`"))
    rows = Vector{Any}[]
    workqueue = [(m, "")]
    parents = Module[m]
    while !isempty(workqueue)
        m2, prep = popfirst!(workqueue)
        for v in names(m2; all, imported)
            if !isdefined(m2, v) || !occursin(pattern, string(v))
                continue
            end
            value = getfield(m2, v)
            isbuiltin = value === Base || value === Main || value === Core
            if recursive && !isbuiltin && isa(value, Module) && value !== m2 && nameof(value) === v && value ∉ parents
                push!(parents, value)
                push!(workqueue, (value, "$v."))
            end
            ssize_str, ssize = if isbuiltin
                ("", typemax(Int))
            else
                ss = summarysize(value)
                (format_bytes(ss), ss)
            end
            if ssize >= minsize
                push!(rows, Any[string(prep, v), ssize_str, summary(value), ssize])
            end
        end
    end
    let (col, rev) = if sortby === :name
            1, false
        elseif sortby === :size
            4, true
        elseif sortby === :summary
            3, false
        else
            @assert "unreachable"
        end
        sort!(rows; by=r -> r[col], rev)
    end
    pushfirst!(rows, Any["name", "size", "summary"])

    return Markdown.MD(Any[Markdown.Table(map(r -> r[1:3], rows), Symbol[:l, :r, :l])])
end

"""
    malloc_trim_if_glibc()

Check if the underlying system uses the GNU C Library (GLIBC) and, if so, calls `malloc_trim(0)`
to attempt to release memory back to the system. This function is primarily useful on Linux systems
where GLIBC is the standard C library.

On non-Linux systems, or Linux systems not using GLIBC, the function does nothing

## Notes

  - This method uses a heuristic by checking for the presence of "glibc" in `/proc/version`.
    While this is indicative of GLIBC's presence on many typical systems, it's not a guaranteed check.
  - The actual `malloc_trim` function is specific to the glibc implementation of the C standard library.
"""
function malloc_trim_if_glibc()
    # Check if on Linux
    if Sys.islinux()
        # Check for glibc by reading /proc/version (this is specific to Linux)
        try
            version_info = read("/proc/version", String)
            if occursin("glibc", version_info)
                ccall(:malloc_trim, Cvoid, (Cint,), 0)
                #println("malloc_trim called.")
            else
                #println("Not using glibc.")
            end
        catch e
            if isa(e, InterruptException)
                rethrow(e)
            end
            #println("Error reading /proc/version: ", e)
        end
    else
        #println("Not on a Linux system.")
    end
end

"""
    extract_dds_to_dataframe(dds::Vector{IMAS.dd{Float64}}, xtract=IMAS.ExtractFunctionsLibrary)

Extracts scalars quantities from ExtractFunctionsLibrary in parallel and return the dataframe
"""
function extract_dds_to_dataframe(dds::Vector{IMAS.dd{Float64}}, xtract=IMAS.ExtractFunctionsLibrary)
    extr = Dict(extract(dds[1], xtract))
    df = DataFrames.DataFrame(extr)
    for k in 2:length(dds)
        push!(df, df[1, :])
    end
    p = ProgressMeter.Progress(length(dds); showspeed=true)
    Threads.@threads for k in eachindex(dds)
        try
            tmp = Dict(extract(dds[k], xtract))
            df[k, :] = tmp
        catch e
            if isa(e, InterruptException)
                rethrow(e)
            end
            continue
        end
        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)
    return df
end

"""
    install_fusebot(folder::String=dirname(readchomp(`which juliaup`)))

This function installs the `fusebot` executable in a given folder,
by default in the directory where the juliaup executable is located.
"""
function install_fusebot(folder::String=dirname(readchomp(`which juliaup`)))
    fusebot_path = joinpath(dirname(dirname(pathof(FUSE))), "fusebot")
    target_path = joinpath(folder, "fusebot")
    ptp_target_path = joinpath(folder, "ptp")

    if !isfile(fusebot_path)
        error("The `fusebot` executable does not exist in the FUSE directory!?")
    end

    cp(fusebot_path, target_path; force=true)

    if folder == dirname(readchomp(`which juliaup`))
        println("`fusebot` has been successfully installed in the Julia executable directory: $folder")
    else
        println("`fusebot` has been successfully installed in folder: $folder")
    end

    if isfile(ptp_target_path)
        rm(ptp_target_path)
        println("Old `ptp` has been successfully removed from folder: $folder")
    end
end
