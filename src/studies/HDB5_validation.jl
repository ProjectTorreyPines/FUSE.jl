import Random
import Dates
import CSV
import Distributed: pmap

#= ====================== =#
#  StudyHDB5Validation     #
#= ====================== =#

"""
    study_parameters(::Val{:HDB5Validation})

Runs FUSE stationary plasma simulations against H-mode confinement database (HDB5)
entries and compares predicted τ_energy against experimental values.

See https://osf.io/593q6/ for database description.
"""
function study_parameters(::Val{:HDB5Validation})
    return FUSEparameters__ParametersStudyHDB5Validation{Real}()
end

Base.@kwdef mutable struct FUSEparameters__ParametersStudyHDB5Validation{T<:Real} <: ParametersStudy{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StudyHDB5Validation
    server::Switch{String} = study_common_parameters(; server="localhost")
    n_workers::Entry{Int} = study_common_parameters(; n_workers=2)
    release_workers_after_run::Entry{Bool} = study_common_parameters(; release_workers_after_run=true)
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the validation runs into")
    tokamak::Entry{Symbol} = Entry{Symbol}("-", "Tokamak to validate against (:all for all tokamaks)"; default=:all)
    n_samples_per_tokamak::Entry{Int} = Entry{Int}("-", "Number of random samples per tokamak (0 for all cases)"; default=10)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat0quench, :sat1, :sat1geo, :sat2, :sat3], "-", "Saturation rule for TGLF"; default=:sat0quench)
    tglfnn_model::Entry{String} = Entry{String}("-", "TGLFNN model name, or \"auto\" to auto-select per tokamak via model_selector"; default="auto")
    n_model_selector_samples::Entry{Int} = Entry{Int}("-", "Number of random cases per tokamak used for model_selector (0 for all available cases)"; default=5)
end

mutable struct StudyHDB5Validation{T<:Real} <: AbstractStudy
    sty::OverrideParameters{T,FUSEparameters__ParametersStudyHDB5Validation{T}}
    act::Union{ParametersAllActors,Missing}
    run_df::Union{DataFrame,Missing}
    fail_df::Union{DataFrame,Missing}
end

function StudyHDB5Validation(sty::ParametersStudy; kw...)
    sty = OverrideParameters(sty; kw...)
    study = StudyHDB5Validation(sty, missing, missing, missing)
    parallel_environment(sty.server, sty.n_workers)
    return study
end

function StudyHDB5Validation(sty::ParametersStudy, act::ParametersAllActors; kw...)
    sty = OverrideParameters(sty; kw...)
    study = StudyHDB5Validation(sty, act, missing, missing)
    parallel_environment(sty.server, sty.n_workers)
    return study
end

"""
    _run(study::StudyHDB5Validation)

Runs the HDB5 validation study: loads HDB5 cases, runs ActorStationaryPlasma on each,
and compares predicted τ_energy with experimental values.
"""
function _run(study::StudyHDB5Validation)
    sty = study.sty

    @assert (sty.n_workers == 0 || sty.n_workers == length(Distributed.workers())) "The number of workers = $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"

    # Load HDB5 database
    run_df = load_hdb5(sty.tokamak)

    # Sample cases if requested
    if sty.n_samples_per_tokamak > 0
        tok_list = unique(run_df[:, "TOK"])
        tok_count = [nrow(run_df[run_df.TOK.==tok, :]) for tok in tok_list]
        run_df = run_df[Random.shuffle(1:nrow(run_df)), :]
        run_df = vcat([run_df[run_df.TOK.==tok, :][1:min(sty.n_samples_per_tokamak, tok_count[idx]), :] for (idx, tok) in enumerate(tok_list)]...)
    end

    # Resolve TGLFNN model per unique TOK (lookup table first, then model_selector for unknown TOKs)
    tok_models = Dict{String,String}()
    unique_toks = unique(run_df.TOK)
    if sty.tglfnn_model == "auto"
        toks_to_search = String[]
        for tok in unique_toks
            lookup_key = (tok, sty.sat_rule)
            if haskey(KNOWN_TGLFNN_MODELS, lookup_key)
                tok_models[tok] = KNOWN_TGLFNN_MODELS[lookup_key]
                @info "Using known TGLFNN model" tok model = tok_models[tok] sat_rule = sty.sat_rule
            else
                push!(toks_to_search, tok)
            end
        end

        if !isempty(toks_to_search)
            @info "Running model_selector for unknown TOKs" toks = toks_to_search sat_rule = sty.sat_rule n_samples = sty.n_model_selector_samples
            rho_grid_sel = range(0.1, 0.9, 9)
            for tok in toks_to_search
                try
                    tok_data = load_hdb5(Symbol(tok))
                    n_avail = nrow(tok_data)
                    n_sel = sty.n_model_selector_samples > 0 ? min(sty.n_model_selector_samples, n_avail) : n_avail
                    sample_indices = Random.shuffle(1:n_avail)[1:n_sel]

                    all_input_tglfs = TurbulentTransport.InputTGLF{Float64}[]
                    for idx in sample_indices
                        try
                            ini_tmp, act_tmp = case_parameters(tok_data[idx, :]; verbose=false)
                            act_tmp.ActorTGLF.tglfnn_model = "sat1_em_d3d"
                            dd_tmp = IMAS.dd()
                            init(dd_tmp, ini_tmp, act_tmp)
                            rho_vec = collect(rho_grid_sel)
                            input_tglfs = TurbulentTransport._unwrap_input_tglfs(
                                TurbulentTransport.InputTGLF(dd_tmp, rho_vec, sty.sat_rule, true, false))
                            append!(all_input_tglfs, input_tglfs)
                        catch e
                            @debug "Skipping case during model selection" tok idx exception = e
                        end
                    end

                    if !isempty(all_input_tglfs)
                        ms = TurbulentTransport.model_selector(all_input_tglfs;
                            filter_sat_rule=sty.sat_rule, max_models=1, ground_truth=true, verbose=false)
                        tok_models[tok] = ms.rankings[1].top_models[1]
                        @info "Auto-selected TGLFNN model" tok model = tok_models[tok] sat_rule = sty.sat_rule n_cases_used = length(sample_indices) n_inputs = length(all_input_tglfs)
                    else
                        tok_models[tok] = "sat1_em_d3d"
                        @warn "No cases initialized successfully, using fallback" tok fallback = tok_models[tok]
                    end
                catch e
                    tok_models[tok] = "sat1_em_d3d"
                    @warn "model_selector failed, using fallback" tok fallback = tok_models[tok] exception = (e, catch_backtrace())
                end
            end
        end
    else
        for tok in unique_toks
            tok_models[tok] = sty.tglfnn_model
        end
    end

    n_cases = nrow(run_df)

    # Add output columns
    run_df[:, "TAUTH_fuse"] = zeros(n_cases)
    run_df[:, "T0_fuse"] = zeros(n_cases)
    run_df[:, "tglfnn_model"] = [tok_models[run_df[k, :TOK]] for k in 1:n_cases]
    run_df[:, "error_message"] = fill("", n_cases)

    sat_rule = sty.sat_rule

    # Run in parallel
    @info "Running HDB5 validation" n_cases workers = sty.n_workers server = sty.server sat_rule
    data_rows = ProgressMeter.@showprogress pmap(
        k -> _run_HDB5_case(run_df[k, :], study.act; sat_rule, tglfnn_model=tok_models[run_df[k, :TOK]]),
        1:n_cases
    )
    for k in 1:n_cases
        run_df[k, :] = data_rows[k]
    end

    # Split into success/failure
    study.fail_df = filter(:TAUTH_fuse => isnan, run_df)
    study.run_df = filter(:TAUTH_fuse => !isnan, run_df)

    n_fail = nrow(study.fail_df)
    n_total = n_fail + nrow(study.run_df)
    @info "Validation complete" failed = n_fail total = n_total
    if nrow(study.run_df) > 0
        MRE = round(100 * _mean_relative_error(study.run_df[:, :TAUTH], study.run_df[:, :TAUTH_fuse]); digits=2)
        @info "Mean Relative Error" MRE_percent = MRE
    end

    # Save results
    if !ismissing(getproperty(sty, :save_folder, missing)) && !isempty(sty.save_folder)
        if !isdir(sty.save_folder)
            mkpath(sty.save_folder)
        end
        timestamp = Dates.format(Dates.now(), "yyyy-mm-ddTHH-MM-SS")
        CSV.write(joinpath(sty.save_folder, "HDB5_validation_$(timestamp).csv"), study.run_df)
        if nrow(study.fail_df) > 0
            CSV.write(joinpath(sty.save_folder, "HDB5_validation_$(timestamp)_failed.csv"), study.fail_df)
        end
    end

    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"
    end

    return study
end

"""
    _run_HDB5_case(data_row::DataFrames.DataFrameRow, study_act::Union{ParametersAllActors,Missing}; sat_rule::Symbol=:sat0quench, tglfnn_model::String="sat1_em_d3d")

Run a single HDB5 validation case: init from data_row, run ActorStationaryPlasma,
and extract predicted τ_energy and T0.
"""
function _run_HDB5_case(data_row::DataFrames.DataFrameRow, study_act::Union{ParametersAllActors,Missing};
                        sat_rule::Symbol=:sat0quench, tglfnn_model::String="sat1_em_d3d")
    try
        ini, case_act = case_parameters(data_row; verbose=false)
        act = ismissing(study_act) ? case_act : study_act
        act.ActorTGLF.sat_rule = sat_rule
        act.ActorTGLF.tglfnn_model = tglfnn_model
        dd = IMAS.dd()
        actor_logging(dd, false)
        init!(dd, ini, act)
        ActorStationaryPlasma(dd, act)
        data_row[:TAUTH_fuse] = @ddtime(dd.summary.global_quantities.tau_energy.value)
        data_row[:T0_fuse] = dd.core_profiles.profiles_1d[].electrons.temperature[1]
        data_row[:error_message] = ""
    catch e
        data_row[:TAUTH_fuse] = NaN
        data_row[:T0_fuse] = NaN
        data_row[:error_message] = "$e"
        if isa(e, InterruptException)
            rethrow(e)
        end
    end
    return data_row
end

"""
    load_HDB5_validation(filename::AbstractString)

Load a previously saved HDB5 validation run from a CSV file and return a `StudyHDB5Validation`.
The failed-cases CSV (if present) is inferred from the filename by appending `_failed`.

Example:
```julia
study = FUSE.load_HDB5_validation("/tmp/hdb5/HDB5_validation_2026-01-01T12-00-00.csv")
plot(study)
```
"""
function load_HDB5_validation(filename::AbstractString)
    sty = OverrideParameters(FUSEparameters__ParametersStudyHDB5Validation{Real}())
    study = StudyHDB5Validation(sty, missing, missing, missing)
    study.run_df = CSV.read(filename, DataFrame)
    if "error_message" in names(study.run_df)
        study.run_df[!, :error_message] = coalesce.(study.run_df[!, :error_message], "")
    end
    failed_filename = replace(filename, ".csv" => "_failed.csv")
    if isfile(failed_filename)
        study.fail_df = CSV.read(failed_filename, DataFrame)
        if "error_message" in names(study.fail_df)
            study.fail_df[!, :error_message] = coalesce.(study.fail_df[!, :error_message], "")
        end
    end
    return study
end

function _R_squared(x, y)
    return 1 - sum((x .- y) .^ 2) / sum((x .- sum(x) / length(x)) .^ 2)
end

function _mean_relative_error(x, y)
    return sum(abs, (y .- x) ./ x) / length(x)
end

@recipe function plot_HDB5_validation(study::StudyHDB5Validation)
    df = study.run_df
    if ismissing(df) || nrow(df) == 0
        @warn "No successful results to plot"
        return
    end

    df = df[DataFrames.completecases(df), :]
    df = filter(row -> row[:TAUTH_fuse] > 0.0 && isfinite(row[:TAUTH_fuse]), df)

    x = df[:, :TAUTH]
    y = df[:, :TAUTH_fuse]

    MRE = round(100 * _mean_relative_error(x, y); digits=2)
    R² = round(_R_squared(x, y); digits=2)

    x_ylim = [2.5e-4, 1e1]
    off = sqrt(MRE / 100)

    all_toks = unique(df.TOK)
    colors = Dict("JET" => :red, "D3D" => :blue, "NSTX" => :purple, "JT60U" => :black,
                   "ASDEX" => :darkgreen, "AUG" => :green4, "CMOD" => :orange)

    aspect_ratio --> :equal
    xscale --> :log10
    yscale --> :log10
    xlim --> x_ylim
    ylim --> x_ylim
    legend --> :topleft
    xlabel --> "Experimental τₑ [s]"
    ylabel --> "FUSE τₑ [s]"

    # MRE band lines
    @series begin
        seriestype := :line
        linestyle := :dash
        color := :black
        label := "±$(MRE)% MRE ($(nrow(df)) cases)"
        [x_ylim[1], x_ylim[2]], [off * x_ylim[1], off * x_ylim[2]]
    end

    @series begin
        seriestype := :line
        linestyle := :dash
        color := :black
        primary := false
        [off * x_ylim[1], off * x_ylim[2]], [x_ylim[1], x_ylim[2]]
    end

    # 1:1 line
    @series begin
        seriestype := :line
        color := :black
        label := nothing
        [x_ylim[1], x_ylim[2]], [x_ylim[1], x_ylim[2]]
    end

    has_model_col = "tglfnn_model" in names(df)

    # Scatter per tokamak
    for tok in all_toks
        index = df.TOK .== tok
        c = get(colors, tok, :gray)
        tok_label = if has_model_col
            models = unique(df[index, :tglfnn_model])
            "$tok ($(join(models, ", ")))"
        else
            "$tok"
        end
        @series begin
            seriestype := :scatter
            color := c
            label := tok_label
            markersize := 2.0
            markerstrokewidth := 0.0
            legend_foreground_color := :transparent
            legend_background_color := :transparent
            df[index, :TAUTH], df[index, :TAUTH_fuse]
        end
    end
end
