import Random
import DataFrames
import CSV
import ProgressMeter

"""
    simple_equilibrium_transport_workflow(dd::IMAS.dd,
                                ini::InitParameters;
                                act::ActorParameters;
                                save_directory::String,
                                do_plot :: Bool)

Initializes and runs simple equilibrium, core_sources and transport actors and stores the resulting dd in <save_directory>
"""
function simple_equilibrium_transport_workflow(dd::IMAS.dd, ini::InitParameters, act::ActorParameters; save_directory::String="", do_plot::Bool=false, warn_nn_train_bounds=true, transport_model=:tglfnn, verbose=false)
    FUSE.init_equilibrium(dd, ini, act) # already solves the equilibrium once
    FUSE.init_core_profiles(dd, ini, act)
    FUSE.init_core_sources(dd, ini, act)

    # Set j_ohmic to steady state
    IMAS.j_ohmic_steady_state!(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])

    # run transport actor
    FUSE.TauennActor(dd, act; transport_model=transport_model, warn_nn_train_bounds, verbose=verbose)

    # Set beta_normal from equilbrium to the kinetic beta_n
    if !isempty(dd.core_profiles.profiles_1d)
        dd.equilibrium.time_slice[].global_quantities.beta_normal = @ddtime(dd.summary.global_quantities.beta_tor_thermal_norm.value)
    end

    # run equilibrium actor with the updated beta
    FUSE.SolovevActor(dd, act)

    if do_plot
        display(plot(dd.equilibrium, psi_levels_out=[], label=ini.general.casename))
        display(plot(dd.core_profiles))
        display(plot(dd.core_sources))
    end

    if !isempty(save_directory)
        IMAS.imas2json(dd, joinpath(save_directory, "$(ini.general.casename).json"))
    end
    return dd
end

"""
    function transport_validation_workflow(
        n_samples_per_tokamak::Union{Integer,Symbol}, # setting this to :all runs the whole database
        save_directory::String,
        show_dd_plots::Bool,
        plot_database::Bool)

Runs n_samples of the HDB5 database (https://osf.io/593q6) and stores results in save_directory
"""
function transport_validation_workflow(;
    tokamak::Union{String,Symbol}=:all,
    n_samples_per_tokamak::Union{Integer,Symbol}=10,
    save_directory::String="",
    show_dd_plots=false,
    plot_database=true,
    verbose=false)

    # load HDB5 database
    run_df = load_hdb5(tokamak)

    # pick cases at random
    if n_samples_per_tokamak !== :all
        tok_list = unique(run_df[:, "TOK"])
        run_df = DataFrames.reduce(
            vcat, [run_df[run_df.TOK.==tok, :][Random.shuffle(1:DataFrames.nrow(run_df[run_df.TOK.==tok, :]))[1:minimum([n_samples_per_tokamak, length(run_df[run_df.TOK.==tok, :][:, "TOK"])])], :] for tok in tok_list]
        )
    end
    
    # Run simple_equilibrium_transport_workflow on each of the selected cases
    tau_FUSE = zeros(length(run_df[:,"TOK"]))
    tbl = DataFrames.Tables.rowtable(run_df)
    failed_runs_ids = Int[]
    p = ProgressMeter.Progress(length(DataFrames.Tables.rows(tbl)); showspeed=true)
    Base.Threads.@threads for idx in 1:length(DataFrames.Tables.rows(tbl))
        try
            dd = IMAS.dd()
            ini = InitParameters(run_df[idx, :])
            act = ActorParameters()
            simple_equilibrium_transport_workflow(dd, ini, act; save_directory, do_plot=show_dd_plots, warn_nn_train_bounds=false)
            tau_FUSE[idx] = @ddtime(dd.summary.global_quantities.tau_energy.value)
            if verbose
                display(println("τ_fuse = $(@ddtime(dd.summary.global_quantities.tau_energy.value)) τ_hdb5 = $(run_df[idx, :TAUTH])"))
            end
        catch e
            push!(failed_runs_ids, idx)
            if verbose
                display(@show e)
            end
        end
        ProgressMeter.next!(p)
    end
    println("Failed runs: $(length(failed_runs_ids)) out of $(length(run_df[:,"TOK"]))")
    run_df[:,"TAUTH_fuse"] = tau_FUSE

    failed_df = run_df[failed_runs_ids, :]

    # save all input data as well as predicted tau to CSV file
    if !isempty(save_directory)
        CSV.write(joinpath(save_directory, "dataframe.csv"), run_df)
        CSV.write(joinpath(save_directory, "failed_runs_dataframe.csv"), failed_df)
    end

    if plot_database
        plot_x_y_regression(run_df, "TAUTH")
    end

    return run_df, failed_df
end

function R_squared(x, y)
    return 1 - sum((x .- y) .^ 2) / sum((x .- sum(x) / length(x)) .^ 2)
end

function mean_relative_error(x, y)
    return sum(abs.((y .- x) ./ x)) / length(x)
end

"""
    plot_x_y_regression(dataframe::DataFrames.DataFrame, name::Union{String,Symbol}="TAUTH")

Plot regression in log-form on x_name and y_name in the dataframe
"""
function plot_x_y_regression(dataframe::DataFrames.DataFrame, name::Union{String,Symbol}="TAUTH")
    x_name = name
    y_name = "$(name)_fuse"
    if x_name == "TAUTH"
        x_ylim = [5e-3, 1e1]
    else
        x_ylim = [minimum(abs.(dataframe[:, x_name])) / 1e1, maximum(dataframe[:, x_name]) * 1e1]
    end
    dataframe = dataframe[DataFrames.completecases(dataframe), :]
    R² = round(R_squared(dataframe[:, x_name], dataframe[:, y_name]), digits=2)
    MRE = round(100 * mean_relative_error(dataframe[:, x_name], dataframe[:, y_name]), digits=2)
    p = plot(dataframe[:, x_name], dataframe[:, y_name], seriestype=:scatter, xaxis=:log, yaxis=:log, ylim=x_ylim, xlim=x_ylim, xlabel=x_name, ylabel=y_name, label="mean_relative_error = $MRE %")
    plot!([0.5 * x_ylim[1], 0.5 * x_ylim[2]], [2 * x_ylim[1], 2 * x_ylim[2]], linestyle=:dash, label="+50%")
    plot!([2 * x_ylim[1], 2 * x_ylim[2]], [0.5 * x_ylim[1], 0.5 * x_ylim[2]], linestyle=:dash, label="-50%", legend=:topleft)
    display(plot!([x_ylim[1], x_ylim[2]], [x_ylim[1], x_ylim[2]], label=nothing))
    println("R² = $(R²), mean_relative_error = $MRE)")
    return p
end

function plot_x_y_regression(filename::String, name::Union{String,Symbol}="TAUTH")
    dataframe = CSV.read(filename, DataFrames.DataFrame)
    plot_x_y_regression(dataframe, name)
end