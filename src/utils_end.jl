import Weave

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
        DD::Vector{<:Union{AbstractString,IMAS.dd}},
        xtract::AbstractDict{Symbol,IMAS.ExtractFunction}=IMAS.ExtractFunctionsLibrary;
        filter_invalid::Symbol=:none)::DataFrames.DataFrame

Extract data from multiple folders or `dd`s and return results in DataFrame format.

Filtering can by done by `:cols` that have all NaNs, `:rows` that have any NaN, both with `:all`, or `:none`.
"""
function IMAS.extract(
    DD::Vector{<:Union{AbstractString,IMAS.dd}},
    xtract::AbstractDict{Symbol,IMAS.ExtractFunction}=IMAS.ExtractFunctionsLibrary;
    filter_invalid::Symbol=:none)::DataFrames.DataFrame

    # test filter_invalid
    @assert filter_invalid in [:none, :cols, :rows, :all] "filter_invalid can only be one of [:none, :cols, :rows, :all]"

    # load the data
    p = ProgressMeter.Progress(length(DD); showspeed=true)
    all_extracts = Dict()
    Threads.@threads for k in eachindex(DD)
        try
            tmp = Dict(extract(DD[k], xtract))
            tmp[:dir] = DD[k]
            all_extracts[k] = tmp
        catch
            continue
        end
        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)

    # to dataframe
    df = DataFrames.DataFrame(all_extracts[1])
    for k in 1:length(DD)
        if k in keys(all_extracts)
            push!(df, pop!(all_extracts, k))
        end
    end

    # filter
    if filter_invalid ∈ [:cols, :all]
        # drop columns that have all NaNs
        isnan_nostring(x::Any) = (typeof(x) <: Number) ? isnan(x) : false
        visnan(x::Vector) = isnan_nostring.(x)
        df = df[:, .!all.(visnan.(eachcol(df)))]
    end
    if filter_invalid ∈ [:rows, :all]
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
        dd::IMAS.dd,
        ini::ParametersAllInits,
        act::ParametersAllActors,
        savedir::AbstractString;
        freeze::Bool=true,
        format::Symbol=:json)

Save FUSE (dd, ini, act) to dd.json/h5, ini.json, and act.json files
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

Save FUSE (dd, ini, act) to dd.json/h5, ini.json, and act.json files and exception stacktrace to "error.txt"
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
    load(savedir::AbstractString; load_dd::Bool=true, load_ini::Bool=true, load_act::Bool=true, skip_on_error::Bool=false)

Read (dd, ini, act) to dd.json/h5, ini.json, and act.json files.

Returns `missing` for files are not there or if `error.txt` file exists in the folder.
"""
function load(savedir::AbstractString; load_dd::Bool=true, load_ini::Bool=true, load_act::Bool=true, skip_on_error::Bool=false)
    if isfile(joinpath(savedir, "error.txt")) && skip_on_error
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

NOTE: `section` is used internally to produce digest PDFs
"""
function digest(
    dd::IMAS.dd;
    terminal_width::Int=136,
    line_char::Char='─',
    section::Int=0)

    #NOTE: this function is defined in FUSE and not IMAS because it uses Plots.jl and not BaseRecipies.jl
    #      also it references ini and act

    sec = 1
    if section ∈ [0, sec]
        IMAS.print_tiled(extract(dd); terminal_width, line_char)
    end

    # equilibrium with build and PFs
    sec += 1
    if !isempty(dd.equilibrium.time_slice) && section ∈ [0, sec]
        println('\u200B')
        p = plot(dd.equilibrium, legend=false)
        if !isempty(dd.build.layer)
            plot!(p[1], dd.build, legend=false)
        end
        if !isempty(dd.pf_active.coil)
            plot!(p[1], dd.pf_active, legend=false, colorbar=false)
        end
        display(p)
    end

    # build layers
    sec += 1
    if !isempty(dd.build.layer) && section ∈ [0, sec]
        println('\u200B')
        display(dd.build.layer)
    end

    # core profiles
    sec += 1
    if !isempty(dd.core_profiles.profiles_1d) && section ∈ [0, sec]
        println('\u200B')
        display(plot(dd.core_profiles, only=1))
    end
    sec += 1
    if !isempty(dd.core_profiles.profiles_1d) && section ∈ [0, sec]
        println('\u200B')
        display(plot(dd.core_profiles, only=2))
    end
    sec += 1
    if !isempty(dd.core_profiles.profiles_1d) && section ∈ [0, sec]
        println('\u200B')
        display(plot(dd.core_profiles, only=3))
    end

    # core sources
    sec += 1
    if !isempty(dd.core_sources.source) && section ∈ [0, sec]
        println('\u200B')
        display(plot(dd.core_sources, only=1))
    end
    sec += 1
    if !isempty(dd.core_sources.source) && section ∈ [0, sec]
        println('\u200B')
        display(plot(dd.core_sources, only=2))
    end
    sec += 1
    if !isempty(dd.core_sources.source) && section ∈ [0, sec]
        println('\u200B')
        display(plot(dd.core_sources, only=3))
    end
    sec += 1
    if !isempty(dd.core_sources.source) && section ∈ [0, sec]
        println('\u200B')
        display(plot(dd.core_sources, only=4))
    end

    # neutron wall loading
    sec += 1
    if !isempty(dd.neutronics.time_slice) && section ∈ [0, sec]
        println('\u200B')
        xlim = extrema(dd.neutronics.first_wall.r)
        xlim = (xlim[1] - ((xlim[2] - xlim[1]) / 10.0), xlim[2] + ((xlim[2] - xlim[1]) / 10.0))
        display(plot(dd.neutronics.time_slice[].wall_loading; xlim))
    end

    # center stack stresses
    sec += 1
    if !ismissing(dd.solid_mechanics.center_stack.grid, :r_oh) && section ∈ [0, sec]
        println('\u200B')
        display(plot(dd.solid_mechanics.center_stack.stress))
    end

    # balance of plant (cannot be plotted right now plotting can only be done when running actor and not from data in dd)
    # sec += 1
    # if !missing(dd.balance_of_plant, :Q_plant) && section ∈ [0, sec]
    # println('\u200B')
    #     display(plot(dd.balance_of_plant))
    # end

    # costing
    sec += 1
    if !ismissing(dd.costing.cost_direct_capital, :cost) && (dd.costing.cost_direct_capital.cost != 0) && section ∈ [0, sec]
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
    outfilename = joinpath(pwd(), "$(replace(title," "=>"_")).pdf")
    tmpdir = mktempdir()
    logger = SimpleLogger(stderr, Logging.Warn)
    try
        filename = redirect_stdout(Base.DevNull()) do
            filename = with_logger(logger) do
                Weave.weave(joinpath(@__DIR__, "digest.jmd");
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
        cp(filename, outfilename, force=true)
        return outfilename
    catch e
        println("Generation of $(basename(outfilename)) failed. See directory: $tmpdir")
    else
        rm(tmpdir, recursive=true, force=true)
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
        "plasma aspect ratio changed" => :aspect_ratio_change,
        "Unable to blend the core-pedestal" => :blend_core_ped,
        "Bad expression" => :bad_expression,
        "Exceeded limits" => :exceed_lim_A,
        "Some stability models have breached their limit threshold:" => :exceed_lim_B,
        "TaskFailedException" => :task_exception,
        "Could not trace closed flux surface" => :flux_surfaces_A,
        "Flux surface at ψ=" => :flux_surfaces_B,
        "stainless_steel.yield_strength" => :CS_stresses,
        "DomainError with" => :Solovev,
        "BoundsError: attempt to access" => :flux_surfaces_C)
    merge!(error_messages, extra_error_messages)

    # go through directories
    for dir in dirs
        filename = joinpath([dir, "error.txt"])
        if !isfile(filename)
           continue
        end
        first_line,second_line = open(filename, "r") do f
           (readline(f),readline(f))
        end
        found = false
        for (err,cat) in error_messages
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
            if show_first_line
                println(first_line)
                println(second_line)
                println(dir)
                println()
            end
        end
    end

    if do_plot
        labels = collect(keys(errors))
        v = collect(map(length,values(errors)))
        index = sortperm(v)[end:-1:1]
        display(pie(["$(rpad(string(length(errors[cat])),8))   $(string(cat))" for cat in labels[index]], v[index], legend=:outerright))
    end

    return errors
end
