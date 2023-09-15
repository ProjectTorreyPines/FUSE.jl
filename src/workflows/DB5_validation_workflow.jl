import Random
import DataFrames
import CSV
import Dates
import Distributed
import Distributed: pmap
import ProgressMeter
ProgressMeter.ijulia_behavior(:clear)

"""
    workflow_simple_stationary_plasma(
        ini::ParametersAllInits,
        act::ParametersAllActors;
        save_directory::String="",
        do_plot::Bool=false,
        warn_nn_train_bounds=true,
        transport_model=:tglfnn,
        verbose=false)

Initializes and runs simple equilibrium, core_sources and transport actors and stores the resulting dd in <save_directory>
"""
function workflow_simple_stationary_plasma(
    ini::ParametersAllInits,
    act::ParametersAllActors;
    save_directory::String="",
    do_plot::Bool=false,
    warn_nn_train_bounds=true,
    transport_model=:tglfnn,
    verbose=false)

    dd = IMAS.dd()
    init_equilibrium!(dd, ini, act) # already solves the equilibrium once
    init_core_profiles!(dd, ini, act)
    init_core_sources!(dd, ini, act)

    act.ActorTauenn.warn_nn_train_bounds = warn_nn_train_bounds
    act.ActorTauenn.transport_model = transport_model
    act.ActorTauenn.verbose = verbose
    ActorStationaryPlasma(dd, act; do_plot=do_plot)

    if !isempty(save_directory)
        IMAS.imas2json(dd, joinpath(save_directory, "$(ini.general.casename).json"))
    end
    return dd
end

"""
    workflow_HDB5_validation(;
        tokamak::Union{String,Symbol}=:all,
        n_samples_per_tokamak::Union{Integer,Symbol}=10,
        save_directory::String="",
        show_dd_plots=false,
        plot_database=true,
        verbose=false,
        act=missing)

Runs n_samples of the HDB5 database (https://osf.io/593q6) and stores results in save_directory
"""
function workflow_HDB5_validation(;
    tokamak::Union{String,Symbol}=:all,
    n_samples_per_tokamak::Union{Integer,Symbol}=10,
    save_directory::String="",
    show_dd_plots=false,
    plot_database=true,
    verbose=false,
    act=missing)

    # load HDB5 database
    run_df = load_hdb5(tokamak)

    # pick cases at random
    if n_samples_per_tokamak !== :all
        tok_list = unique(run_df[:, "TOK"])
        run_df = DataFrames.reduce(
            vcat,
            [
                run_df[run_df.TOK.==tok, :][
                    Random.shuffle(1:DataFrames.nrow(run_df[run_df.TOK.==tok, :]))[1:minimum([n_samples_per_tokamak, length(run_df[run_df.TOK.==tok, :][:, "TOK"])])],
                    :
                ] for tok in tok_list
            ]
        )
    end

    n_cases = length(run_df.TOK)

    # Outputs to store
    run_df[:, "TAUTH_fuse"] = zeros(n_cases)
    run_df[:, "T0_fuse"] = zeros(n_cases)
    run_df[:, "error_message"] = ["" for i in 1:n_cases]

    # Run workflow_simple_stationary_plasma on each of the selected case
    data_rows = ProgressMeter.@showprogress pmap(row -> run_HDB5_from_data_row(row, act, verbose, show_dd_plots), [run_df[k, :] for k in 1:n_cases])

    for k in 1:length(data_rows)
        run_df[k, :] = data_rows[k]
    end

    failed_df = filter(:TAUTH_fuse => isnan, run_df)
    run_df = filter(:TAUTH_fuse => !isnan, run_df)

    println("Failed runs: $(length(failed_df.TOK)) out of $(length(run_df.TOK))")
    println("Mean Relative error $(MRE = round(100 * mean_relative_error(run_df[:, :TAUTH], run_df[:, :TAUTH_fuse]), digits=2))%")

    # save all input data as well as predicted tau to CSV file
    if !isempty(save_directory)
        CSV.write(joinpath(save_directory, "dataframe_$(Dates.now()).csv"), run_df)
        CSV.write(joinpath(save_directory, "failed_runs_dataframe_$(Dates.now()).csv"), failed_df)
    end

    if plot_database
        plot_x_y_regression(run_df, "TAUTH")
    end

    return run_df, failed_df
end

function run_HDB5_from_data_row(data_row, act::Union{ParametersActor,Missing}=missing, verbose::Bool=false, do_plot::Bool=false)
    try
        ini, ACT = case_parameters(data_row)
        if ismissing(act)
            act = ACT
        end
        dd = workflow_simple_stationary_plasma(ini, act; warn_nn_train_bounds=verbose, do_plot=do_plot)
        data_row[:TAUTH_fuse] = @ddtime (dd.summary.global_quantities.tau_energy.value)
        data_row[:T0_fuse] = dd.core_profiles.profiles_1d[].electrons.temperature[1]
        data_row[:error_message] = ""
    catch e
        data_row[:TAUTH_fuse] = NaN
        data_row[:T0_fuse] = NaN
        data_row[:error_message] = "$e"
    end
    return data_row
end

"""
    plot_x_y_regression(dataframe::DataFrames.DataFrame, name::Union{String,Symbol}="TAUTH")

Plot regression of `\$name` and `\$(name)_fuse` data stored in a given dataframe
"""
function plot_x_y_regression(dataframe::DataFrames.DataFrame, name::Union{String,Symbol}="TAUTH")

    x_name = name
    y_name = "$(name)_fuse"
    if x_name == "TAUTH"
        x_ylim = [5e-3, 1e1]
    else
        x_ylim = [minimum(abs, dataframe[:, x_name]) / 1e1, maximum(dataframe[:, x_name]) * 1e1]
    end
    dataframe = dataframe[DataFrames.completecases(dataframe), :]
    dataframe = filter(row -> row["$(name)_fuse"] > 0.0, dataframe)

    R² = round(R_squared(dataframe[:, x_name], dataframe[:, y_name]); digits=2)
    MRE = round(100 * mean_relative_error(dataframe[:, x_name], dataframe[:, y_name]); digits=2)
    p = plot(
        dataframe[:, x_name],
        dataframe[:, y_name];
        seriestype=:scatter,
        xaxis=:log,
        yaxis=:log,
        ylim=x_ylim,
        xlim=x_ylim,
        xlabel=x_name,
        ylabel=y_name,
        label="mean_relative_error = $MRE % for N = $(length(dataframe[:, x_name]))"
    )
    plot!([0.5 * x_ylim[1], 0.5 * x_ylim[2]], [2 * x_ylim[1], 2 * x_ylim[2]]; linestyle=:dash, label="+50%")
    plot!([2 * x_ylim[1], 2 * x_ylim[2]], [0.5 * x_ylim[1], 0.5 * x_ylim[2]]; linestyle=:dash, label="-50%", legend=:topleft)
    display(plot!([x_ylim[1], x_ylim[2]], [x_ylim[1], x_ylim[2]]; label=nothing))

    println("R² = $(R²), mean_relative_error = $MRE)")
    return p
end

"""
    function plot_x_y_regression(filename::String, name::Union{String,Symbol}="TAUTH")

Plot regression of `\$name` and `\$(name)_fuse` data stored in a given CSV file
"""
function plot_x_y_regression(filename::String, name::Union{String,Symbol}="TAUTH")
    dataframe = CSV.read(filename, DataFrames.DataFrame)
    return plot_x_y_regression(dataframe, name)
end

function R_squared(x, y)
    return 1 - sum((x .- y) .^ 2) / sum((x .- sum(x) / length(x)) .^ 2)
end

function mean_relative_error(x, y)
    return sum(abs, (y .- x) ./ x) / length(x)
end