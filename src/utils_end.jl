import Weave
using InteractiveUtils: summarysize, format_bytes, Markdown
import DelimitedFiles
import OrderedCollections
import DataFrames
import Dates
import HDF5

# ========== #
# Checkpoint #
# ========== #
"""
Provides handy checkpoint capabilities of `dd`, `ini`, `act`

@checkin and @checkout macros use FUSE.checkpoint = Checkpoint()

    ...init...
    @checkin :init dd ini act # store dd, ini, act

    --------

    @checkout :init dd ini act # restore dd, ini, act
    ...something...
    @checkin :something dd act # store dd, act

    --------

    @checkout :something dd # restore only dd
    ...something_else...
    @checkin :something_else dd # store only dd

    --------

    @checkout :something_else dd # restore only dd
    ...something_else_more...

    ========

    empty!(FUSE.checkpoint) # to start over
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

function Base.iterate(chk::Checkpoint, state=1)
    return iterate(chk.history, state)
end

# Generate the delegated methods (NOTE: these methods do not require deepcopy)
for func in [:empty!, :delete!, :haskey, :pop!, :popfirst!]
    @eval function Base.$func(chk::Checkpoint, args...; kw...)
        return $func(chk.history, args...; kw...)
    end
end

"""
    @checkin :key a b c

Macro to save variables into a Checkpoint under a specific key
"""
macro checkin(key, vars...)
    key = esc(key)

    # Save all the variables in the `vars` list under the provided key using their names
    return quote
        @assert typeof($key) <: Symbol "`@checkin chk :what var1 var2` was deprecated in favor of `@checkin :what var1 var2`"
        d = getfield(checkpoint, :history)
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
    @checkout :key a c

Macro to load variables from a Checkpoint
"""
macro checkout(key, vars...)
    key = esc(key)

    # Restore variables from the checkpoint
    return quote
        @assert typeof($key) <: Symbol "`@checkout chk :what var1 var2` was deprecated in favor of `@checkout :what var1 var2`"
        d = getfield(checkpoint, :history)
        if haskey(d, $key)
            saved_vars = d[$key]
            $(Expr(:block, [:($(esc(v)) = deepcopy(getfield(saved_vars, Symbol($(string(v)))))) for v in vars]...))
        else
            throw(Exception("Checkpoint `$key` does not exist"))
        end
        nothing
    end
end

const checkpoint = Checkpoint()

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
        ProgressMeter.ijulia_behavior(:clear)
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
        freeze::Bool=false,
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
    freeze::Bool=false,
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


function save2hdf(
    savedir::AbstractString,
    parent_group::AbstractString,
    dd::Union{Nothing,IMAS.dd},
    ini::Union{Nothing,ParametersAllInits},
    act::Union{Nothing,ParametersAllActors},
    log_io::IOStream;
    error_info::Any=nothing,
    timer::Bool=true,
    varinfo::Bool=false,
    freeze::Bool=false,
    overwrite_groups::Bool=false,
    verbose::Bool=false,
    kw...
)

    savedir = abspath(savedir)
    if !isdir(savedir)
        mkdir(savedir)
    end

    parent_group = IMAS.norm_hdf5_path(parent_group)

    h5_filename = joinpath(savedir, "pid$(getpid())_output.h5")

    function check_and_create_group(fid::HDF5.File, target_group::AbstractString)
        if haskey(fid, target_group)
            if target_group == "/"
                gparent = fid
            else
                if !overwrite_groups
                    error("Target group '$target_group' already exists in file '$(fid.filename)'. " *
                          "\n       Set `overwrite_groups`=true to replace the existing group.")
                else
                    verbose && @warn "Target group '$target_group' already exists. Overwriting it..."
                    HDF5.delete_object(fid, target_group)
                    gparent = HDF5.create_group(fid, target_group)
                end
            end
        else
            gparent = HDF5.create_group(fid, target_group)
        end
        attr = HDF5.attrs(gparent)
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        return gparent
    end

    function check_and_write(fid::HDF5.File, target_group::AbstractString, data)
        if haskey(fid, target_group)
            if target_group != "/"
                if !overwrite_groups
                    error("Target group '$target_group' already exists in file '$(fid.filename)'. " *
                          "\n       Set `overwrite_groups`=true to replace the existing group.")
                else
                    verbose && @warn "Target group '$target_group' already exists. Overwriting it..."
                    HDF5.delete_object(fid, target_group)
                end
            end
        end
        HDF5.write(fid, target_group, data)
        attr = HDF5.attrs(fid[target_group])
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    end

    mode = isfile(h5_filename) ? "r+" : "w"

    HDF5.h5open(h5_filename, mode) do fid

        attr = HDF5.attrs(fid)
        attr["FUSE_version"] = string(pkgversion(FUSE))
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        attr["original_file_abs_path"] = abspath(fid.filename)
        attr["original_file_rel_path"] = relpath(fid.filename)

        if !haskey(fid, parent_group)
            HDF5.create_group(fid, parent_group)
        end
        attr = HDF5.attrs(fid[parent_group])
        attr["FUSE_version"] = string(pkgversion(FUSE))
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        attr["original_file_abs_path"] = abspath(fid.filename)
        attr["original_file_rel_path"] = relpath(fid.filename)

        # Write error information into the HDF5 file (instead of separate txt file)
        if error_info !== nothing
            error_str = ""
            if typeof(error_info) <: Exception
                io = IOBuffer()
                showerror(io, error_info, catch_backtrace())
                error_str = String(take!(io))
            else
                error_str = string(error_info)
            end
            check_and_write(fid, parent_group * "/error.txt", error_str)
        end

        if ini !== nothing
            gparent = check_and_create_group(fid, parent_group * "/ini.h5")
            SimulationParameters.par2hdf!(ini, gparent)
        end

        if dd !== nothing
            IMAS.imas2hdf(dd, h5_filename; mode="a", freeze, target_group=parent_group * "/dd.h5", overwrite=overwrite_groups, verbose)
        end

        if act !== nothing
            gparent = check_and_create_group(fid, parent_group * "/act.h5")
            SimulationParameters.par2hdf!(act, gparent)
        end

        # save timer output
        if timer
            check_and_write(fid, parent_group * "/timer.txt", string(FUSE.timer))
        end

        # save memory trace
        if parse(Bool, get(ENV, "FUSE_MEMTRACE", "false"))
            memtrace_string = String[]
            for (date, txt, kb) in FUSE.memtrace.data
                push!(memtrace_string, "$date $kb \"$txt\"")
            end
            check_and_write(fid, parent_group * "/memtrace.txt", memtrace_string)
        end

        # save vars usage
        if varinfo
            varinfo_string = string(FUSE.varinfo(FUSE; all=true, imported=true, recursive=true, sortby=:size, minsize=1024))
            check_and_write(fid, parent_group * "/varinfo.txt", varinfo_string)
        end

        # save log
        flush(log_io)
        seekstart(log_io)
        log_str = read(log_io, String)
        if !isempty(log_str)
            check_and_write(fid, parent_group * "/log.txt", log_str)
        end
    end

    return savedir
end

function load_database(filename::AbstractString; kw...)
    @assert HDF5.ishdf5(filename) "\"$filename\" is not the HDF5 format"

    HDF5.h5open(filename, "r") do H5_fid
        df = coalesce.(CSV.read(IOBuffer(H5_fid["/extract.csv"][]), DataFrame), NaN)
        return load_database(filename, df[!, :gparent]; kw...)
    end
end

function load_database(filename::AbstractString, conditions::Function; kw...)
    @assert HDF5.ishdf5(filename) "\"$filename\" is not the HDF5 format"

    HDF5.h5open(filename, "r") do H5_fid
        df = coalesce.(CSV.read(IOBuffer(H5_fid["/extract.csv"][]), DataFrame), NaN)
        parent_groups = filter(conditions, df)[!, :gparent]
        return load_database(filename, parent_groups; kw...)
    end
end

function load_database(filename::AbstractString, parent_group::AbstractString, kw...)
    return load_database(filename, [parent_group]; kw...)
end

function load_database(filename::AbstractString, parent_groups::Vector{<:AbstractString}; pattern::Regex=r"", kw...)
    @assert HDF5.ishdf5(filename) "\"$filename\" is not the HDF5 format"

    parent_groups = IMAS.norm_hdf5_path.(parent_groups)

    H5_fid = HDF5.h5open(filename, "r")

    # Load dataframe (from extract)
    df = coalesce.(CSV.read(IOBuffer(H5_fid["/extract.csv"][]), DataFrame), NaN)
    df = subset(df, :gparent => ByRow(x -> x in parent_groups))

    Nparents = length(parent_groups)

    # Prepare output data
    dds = occursin(pattern, "dd.h5") ? fill(IMAS.dd(), Nparents) : nothing
    inis = occursin(pattern, "ini.h5") ? fill(ParametersInits(), Nparents) : nothing
    acts = occursin(pattern, "act.h5") ? fill(ParametersActors(), Nparents) : nothing
    logs = occursin(pattern, "log.txt") ? fill("", Nparents) : nothing
    timers = occursin(pattern, "timer.txt") ? fill("", Nparents) : nothing
    errors = occursin(pattern, "error.txt") ? fill("", Nparents) : nothing

    for (k, gparent) in pairs(parent_groups)
        filterd_keys = filter(x->occursin(pattern,x), keys(H5_fid[gparent]))
        for key in filterd_keys
            h5path = gparent * "/" * key
            if key == "dd.h5"
                dds[k] = IMAS.hdf2imas(filename, h5path)
            elseif key == "ini.h5"
                inis[k] = SimulationParameters.hdf2par(H5_fid[h5path], ParametersInits())
            elseif key == "act.h5"
                acts[k] = SimulationParameters.hdf2par(H5_fid[h5path], ParametersActors())
            elseif key == "log.txt"
                logs[k] = H5_fid[h5path][]
            elseif key == "timer.txt"
                timers[k] = H5_fid[h5path][]
            elseif key == "error.txt"
                errors[k] = H5_fid[h5path][]
            end
        end
    end

    close(H5_fid)

    return (dds=dds, inis=inis, acts=acts, logs=logs, timers=timers, errors=errors, df=df)
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
    for k in 1:3
        if !isempty(dd.core_profiles.profiles_1d) && section ∈ (0, sec)
            println('\u200B')
            display(plot(dd.core_profiles; only=k))
        end
    end

    # core sources
    for k in 1:5+length(IMAS.list_ions(dd.core_sources, dd.core_profiles; time0=dd.global_time))
        if !isempty(dd.core_sources.source) && section ∈ (0, sec)
            println('\u200B')
            display(plot(dd.core_sources; only=k))
        end
    end

    # core transport
    for k in 1:4+length(IMAS.list_ions(dd.core_transport, dd.core_profiles; time0=dd.global_time))
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
        plot!(p, dd.pf_active; what=:currents, time0, title="PF currents at t=$(time0) s", subplot=1)
        plot!(p, dd.equilibrium; time0, cx=true, subplot=2)
        plot!(p, dd.build; subplot=2, legend=false, equilibrium=false, pf_active=false)
        plot!(p, dd.pf_active; time0, subplot=2, coil_identifiers=true)
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
        filename = joinpath(dir, "error.txt")
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
    if Sys.iswindows()
        pid = getpid()
        # Use PowerShell to get the current process's WorkingSet (memory in bytes)
        cmd = `powershell -Command "(Get-Process -Id $pid).WorkingSet64"`
        mem_bytes_str = readchomp(cmd)
        mem_bytes = parse(Int, mem_bytes_str)
    else
        pid = getpid()
        mem_info = read(`ps -p $pid -o rss=`, String)
        mem_usage_kb = parse(Int, strip(mem_info))
        mem_bytes = mem_usage_kb * 1024
    end
    return mem_bytes::Int
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
    ProgressMeter.ijulia_behavior(:clear)
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
    install_fusebot()

Installs the `fusebot` executable in the directory where the `juliaup` executable is located
"""
function install_fusebot()
    try
        if Sys.iswindows()
            folder = dirname(readchomp(`where juliaup`))
        else
            folder = dirname(readchomp(`which juliaup`))
        end
        return install_fusebot(folder)
    catch e
        error("error locating `juliaup` executable: $(string(e))\nPlease use `FUSE.install_fusebot(folder)` specifying a folder in your \$PATH")
    end
end

"""
    install_fusebot(folder::String)

Installs the `fusebot` executable in a specified folder
"""
function install_fusebot(folder::String)
    fusebot_path = joinpath(dirname(dirname(pathof(FUSE))), "fusebot")
    target_path = joinpath(folder, "fusebot")
    @assert isfile(fusebot_path) "The `fusebot` executable does not exist in the FUSE directory!?"
    cp(fusebot_path, target_path; force=true)
    return println("`fusebot` has been successfully installed: $target_path")
end

"""
    compare_manifests(env1_dir::AbstractString, env2_dir::AbstractString)

This function activates the `Manifest.toml` files for the provided directories and compares their dependencies. It identifies:

  - **Added dependencies**: Packages present in the env2 environment but not in the working environment.
  - **Removed dependencies**: Packages present in the working environment but not in the env2 environment.
  - **Modified dependencies**: Packages that exist in both environments but differ in version.
"""
function compare_manifests(env1_dir::AbstractString, env2_dir::AbstractString)
    # Save the current active environment
    original_env = Base.current_project()

    try
        # Activate the env2 environment and retrieve its dependencies
        Pkg.activate(env2_dir)
        env_env2 = Pkg.dependencies()

        # Activate the working environment and retrieve its dependencies
        Pkg.activate(env1_dir)
        env_env1 = Pkg.dependencies()

        # Compare dependencies
        added = setdiff(keys(env_env2), keys(env_env1))
        removed = setdiff(keys(env_env1), keys(env_env2))
        modified = [uuid for uuid in intersect(keys(env_env1), keys(env_env2)) if env_env1[uuid] != env_env2[uuid]]

        println("Added dependencies: ")
        for uuid in added
            package_pkg = env_env2[uuid]
            package_name = package_pkg.name
            package_version = package_pkg.version
            println("    $package_name: $package_version")
        end
        println()
        println("Removed dependencies:")
        for uuid in removed
            package_pkg = env_env1[uuid]
            package_name = package_pkg.name
            package_version = package_pkg.version
            println("    $package_name: $package_version")
        end
        println()
        println("Modified dependencies:")
        for uuid in modified
            env2_pkg = env_env2[uuid]
            env1_pkg = env_env1[uuid]
            env2_name = env2_pkg.name
            env2_version = env2_pkg.version
            env1_version = env1_pkg.version
            println("    $env2_name: env2=$env2_version, env1=$env1_version")
        end

        return (added=added, removed=removed, modified=modified)
    finally
        # Restore the original environment
        Pkg.activate(original_env)
    end
end
