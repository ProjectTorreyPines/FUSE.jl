import Random
import DataFrames
import CSV
import Statistics

"""
    simple_equilibrium_transport_workflow(dd::IMAS.dd,
                                par::Parameters,
                                save_directory::String,
                                show_dd_plots :: Bool)

Initializes and runs simple equilibrium, core_sources and transport actors and stores the resulting dd in <save_directory>
"""
function simple_equilibrium_transport_workflow(dd, par; save_directory::String=missing, show_dd_plots=false)
    FUSE.init_equilibrium(dd, par)
    FUSE.init_core_profiles(dd, par)
    FUSE.init_core_sources(dd, par)

    # run transport actor
    FUSE.TauennActor(dd, par, transport_model=:tglfnn, verbose=false)

    # Set beta_normal from equilbrium to the kinetic beta_n
    if !isempty(dd.core_profiles.profiles_1d)
        dd.equilibrium.time_slice[].global_quantities.beta_normal = @ddtime(dd.summary.global_quantities.beta_tor_thermal_norm.value)
    end

    # run equilibrium actor
    FUSE.SolovevEquilibriumActor(dd, par)

    # correct equilibrium volume and area
    if !ismissing(par.equilibrium, :volume)
        dd.equilibrium.time_slice[].profiles_1d.volume .*= par.equilibrium.volume / dd.equilibrium.time_slice[].profiles_1d.volume[end]
    end
    if !ismissing(par.equilibrium, :area)
        dd.equilibrium.time_slice[].profiles_1d.area .*= par.equilibrium.area / dd.equilibrium.time_slice[].profiles_1d.area[end]
    end

    if show_dd_plots
        display(plot(dd.equilibrium, psi_levels_out=[], label=par.general.casename))
        display(plot(dd.core_profiles))
        display(plot(dd.core_sources))
    end

    if !isempty(save_directory)
        IMAS.imas2json(dd, joinpath(save_directory, "$(par.general.casename).json"))
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
function transport_validation_workflow(
    n_samples_per_tokamak::Union{Integer,Symbol}=10;
    save_directory::String="",
    show_dd_plots=false,
    plot_database=true)

    # Set up the database to run
    # For description of variables see https://osf.io/593q6/
    run_df = load_dataframe(joinpath(dirname(abspath(@__FILE__)), "..", "..", "sample", "HDB5_compressed.csv"))
    signal_names = ["TOK", "SHOT", "AMIN", "KAPPA", "DELTA", "NEL", "ZEFF", "TAUTH", "RGEO", "BT", "IP", "PNBI", "ENBI", "PICRH", "PECRH", "POHM", "MEFF", "VOL", "AREA", "WTH"]
    # subselect on the signals of interest
    run_df = run_df[:, signal_names]
    # only retain cases for which all signals have data
    run_df = run_df[DataFrames.completecases(run_df), :]
    # some basic filters
    run_df = run_df[(run_df.TOK.!="T10").&(run_df.TOK.!="TDEV").&(run_df.KAPPA.>1.0).&(1.6 .< run_df.MEFF .< 2.2).&(1.1 .< run_df.ZEFF .< 5.9), :]
    # pick cases at random
    if n_samples_per_tokamak !== :all
        tok_list = unique(run_df[:, "TOK"])
        run_df = DataFrames.reduce(
            vcat, [run_df[run_df.TOK.==tok, :][Random.shuffle(1:DataFrames.nrow(run_df[run_df.TOK.==tok, :]))[1:minimum([n_samples_per_tokamak, length(run_df[run_df.TOK.==tok, :][:, "TOK"])])], :] for tok in tok_list]
        )
    end

    # Run simple_equilibrium_transport_workflow on each of the selected cases
    run_df[!, :τ_fuse] = tau_FUSE = Vector{Union{Real,Missing}}(missing, length(run_df[:, "TOK"]))
    tbl = DataFrames.Tables.rowtable(run_df)
    failed_runs_ids = Int[]
    for (idx, data_row) in enumerate(DataFrames.Tables.rows(tbl))
        try
            dd = IMAS.dd()
            par = Parameters(:HDB5; data_row)
            simple_equilibrium_transport_workflow(dd, par; save_directory, show_dd_plots)
            tau_FUSE[idx] = @ddtime(dd.summary.global_quantities.tau_energy.value)
        catch e
            push!(failed_runs_ids, idx)
            @error e
        end
    end
    println("Failed runs: $(length(failed_runs_ids)) out of $(length(run_df[:,"TOK"]))")

    if plot_database
        plot_x_y_regression(run_df, "TAUTH", "τ_fuse")
    end

    # save all input data as well as predicted tau to CSV file
    if !isempty(save_directory)
        CSV.write(joinpath(save_directory, "dataframe.csv"), run_df)
    end

    return run_df
end

function load_dataframe(filename::String)
    return CSV.read(filename, DataFrames.DataFrame)
end

function R_squared(x, y)
    return 1 - sum((x .- y) .^ 2) / sum((x .- Statistics.mean(x)) .^ 2)
end

function mean_relative_error(x, y)
    return sum(abs.((y .- x) ./ x)) / length(x)
end

"""
    plot_x_y_regression(dataframe::DataFrames.DataFrame, x_name::Union{String,Symbol}="TAUTH", y_name::Union{String,Symbol}="τ_fuse")

Plot regression in log-form on x_name and y_name in the dataframe
"""
function plot_x_y_regression(dataframe::DataFrames.DataFrame, x_name::Union{String,Symbol}="TAUTH", y_name::Union{String,Symbol}="τ_fuse")
    if x_name == "TAUTH"
        x_ylim = [5e-3, 1e1]
    else
        x_ylim = [minimum(abs.(dataframe[:, x_name])) / 1e1, maximum(dataframe[:, x_name]) * 1e1]
    end
    DataFrames.dropmissing!(dataframe)
    R² = round(R_squared(dataframe[:, x_name], dataframe[:, y_name]), digits=2)
    MRE = round(100 * mean_relative_error(dataframe[:, x_name], dataframe[:, y_name]), digits=2)
    p = plot(dataframe[:, x_name], dataframe[:, y_name], seriestype=:scatter, xaxis=:log, yaxis=:log, ylim=x_ylim, xlim=x_ylim, xlabel=x_name, ylabel=y_name, label="mean_relative_error = $MRE %")
    plot!([0.5 * x_ylim[1], 0.5 * x_ylim[2]], [2 * x_ylim[1], 2 * x_ylim[2]], linestyle=:dash, label="+50%")
    plot!([2 * x_ylim[1], 2 * x_ylim[2]], [0.5 * x_ylim[1], 0.5 * x_ylim[2]], linestyle=:dash, label="-50%", legend=:topleft)
    display(plot!([x_ylim[1], x_ylim[2]], [x_ylim[1], x_ylim[2]], label=nothing))
    println("R² = $(R²), mean_relative_error = $MRE)")
    return p
end

function plot_x_y_regression(filename::String, x_name::Union{String,Symbol}="TAUTH", y_name::Union{String,Symbol}="τ_fuse")
    dataframe = CSV.read(filename, DataFrames.DataFrame)
    plot_x_y_regression(dataframe, x_name, y_name)
end