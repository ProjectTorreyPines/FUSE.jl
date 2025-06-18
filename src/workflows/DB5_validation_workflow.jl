import Random
import DataFrames
import Serialization
import Dates
import Distributed
import Distributed: pmap

"""
    workflow_HDB5_validation(;
        tokamak::Union{String,Symbol}=:all,
        n_samples_per_tokamak::Union{Integer,Symbol}=10,
        save_directory::String="",
        verbose::Bool=false,
        act::Union{ParametersAllActors,Missing}=missing)

Runs n_samples of the HDB5 database (https://osf.io/593q6) and stores results in save_directory
"""
function workflow_HDB5_validation(;
    tokamak::Union{String,Symbol}=:all,
    n_samples_per_tokamak::Union{Integer,Symbol}=10,
    save_directory::String="",
    verbose::Bool=false,
    act::Union{ParametersAllActors,Missing}=missing)

    # load HDB5 database
    run_df = load_hdb5(tokamak)

    # pick cases at random
    if n_samples_per_tokamak !== :all
        tok_list = unique(run_df[:, "TOK"])
        tok_count = [length(run_df[run_df.TOK.==tok, :].TOK) for tok in tok_list]
        run_df = run_df[DataFrames.shuffle(1:length(run_df.TOK)), :]
        run_df = vcat([run_df[run_df.TOK.==tok, :][1:minimum([n_samples_per_tokamak, tok_count[idx]]), :] for (idx, tok) in enumerate(tok_list)]...)
        display(run_df)
    end

    n_cases = length(run_df.TOK)

    # Outputs to store
    run_df[:, "TAUTH_fuse"] = zeros(n_cases)
    run_df[:, "T0_fuse"] = zeros(n_cases)
    run_df[:, "error_message"] = ["" for i in 1:n_cases]

    # Run workflow_simple_stationary_plasma on each of the selected case
    data_rows = ProgressMeter.@showprogress pmap(data_row -> run_HDB5_from_data_row(data_row, act; verbose), [run_df[k, :] for k in 1:n_cases])
    for k in 1:n_cases
        run_df[k, :] = data_rows[k]
    end

    fail_df = filter(:TAUTH_fuse => isnan, run_df)
    run_df = filter(:TAUTH_fuse => !isnan, run_df)

    println("Failed runs: $(length(fail_df.TOK)) out of $(length(fail_df.TOK) + length(run_df.TOK))")
    println("Mean Relative error $(MRE = round(100 * mean_relative_error(run_df[:, :TAUTH], run_df[:, :TAUTH_fuse]), digits=2))%")

    # save all input data as well as predicted tau to JLS file
    if !isempty(save_directory)
        filename = "HDB5_runs_dataframe_$(Dates.now())"
        Serialization.serialize(joinpath(save_directory, "$(filename).jls"), run_df)
        Serialization.serialize(joinpath(save_directory, "$(filename)_failed.jls"), fail_df)
    end

    return (run_df=run_df, fail_df=fail_df)
end

function run_HDB5_from_data_row(data_row, act::Union{ParametersAllActors,Missing}=missing; verbose::Bool=false)
    try
        ini, ACT = case_parameters(data_row; verbose)
        if ismissing(act)
            act = ACT
        end
        dd = IMAS.dd()
        actor_logging(dd, verbose)
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
        elseif verbose
            @error repr(e)
        end
    end
    return data_row
end

"""
    plot_τ_regression(dataframe::DataFrames.DataFrame; kw...)

Plot HDB5 validation workflow stored in a given dataframe

kw arguments are passed to the plot
"""
function plot_τ_regression(dataframe::DataFrames.DataFrame; kw...)
    x_name = "TAUTH"
    y_name = "TAUTH_fuse"
    xlabel = "Experiments τₑ"
    ylabel = "FUSE τₑ"
    x_ylim = [5e-3, 1e1]

    dataframe = dataframe[DataFrames.completecases(dataframe), :]
    dataframe = filter(row -> row[y_name] > 0.0 && isfinite(row[y_name]), dataframe)

    R² = round(R_squared(dataframe[:, x_name], dataframe[:, y_name]); digits=2)
    MRE = round(100 * mean_relative_error(dataframe[:, x_name], dataframe[:, y_name]); digits=2)

    p = plot(; palette=:tab10, aspect_ratio=:equal)
    plot!(
        dataframe[:, x_name],
        dataframe[:, y_name];
        seriestype=:scatter,
        group=dataframe.TOK,
        xscale=:log10,
        yscale=:log10,
        ylim=x_ylim,
        xlim=x_ylim,
        xlabel=xlabel,
        ylabel=ylabel,
        title="$(length(dataframe[:, x_name])) cases",
        kw...
    )

    off = sqrt(MRE / 100)
    plot!([x_ylim[1], x_ylim[2]], [off * x_ylim[1], off * x_ylim[2]]; linestyle=:dash, label="±$(MRE)%", legend=:topleft, color=:black)
    plot!([off * x_ylim[1], off * x_ylim[2]], [x_ylim[1], x_ylim[2]]; linestyle=:dash, color=:black, primary=false)
    plot!([x_ylim[1], x_ylim[2]], [x_ylim[1], x_ylim[2]]; label=nothing, color=:black)

    println("R² = $(R²), mean_relative_error = $MRE)")
    return p
end

"""
    load_hdb5(filename::String)

Load h5db from as JLS file
"""
function load_hdb5(filename::String)
    return Serialization.deserialize(filename)
end

"""
    plot_τ_regression(filename::String)

Plot τe regression from data stored in a given JLS file
"""
function plot_τ_regression(filename::String; kw...)
    dataframe = load_hdb5(filename)
    return plot_τ_regression(dataframe; kw...)
end

function R_squared(x, y)
    return 1 - sum((x .- y) .^ 2) / sum((x .- sum(x) / length(x)) .^ 2)
end

function mean_relative_error(x, y)
    return sum(abs, (y .- x) ./ x) / length(x)
end