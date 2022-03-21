import Random

include("load_plot_dataframe.jl")

function run_equilibrium_for_database(dd, par)
    # Set beta_normal from equilbrium to the kinetic beta_n
    if !isempty(dd.core_profiles.profiles_1d)
        dd.equilibrium.time_slice[].global_quantities.beta_normal = @ddtime(dd.summary.global_quantities.beta_tor_thermal_norm.value)
    end
    eqactor = FUSE.SolovevEquilibriumActor(dd.equilibrium)
    FUSE.step(eqactor, verbose = false)
    FUSE.finalize(eqactor)

    # correct equilibrium volume and area
    eqt1d = dd.equilibrium.time_slice[].profiles_1d
    if !ismissing(par.equilibrium, :volume)
        eqt1d.volume .*= par.equilibrium.volume / eqt1d.volume[end]
    end
    if !ismissing(par.equilibrium, :area)
        eqt1d.area .*= par.equilibrium.area / eqt1d.area[end]
    end
end

"""
    function validation_workflow(
        n_samples::Union{Real,Symbol},
        save_directory::String,
        show_equilibrium::Bool,
        plot_database::Bool)

Runs n_samples of the HDB5 database and stores results in save_directory
"""
function simple_equilibrium_transport_workflow(dd, par; save_directory::String = missing, show_dd_plots = false)
    FUSE.init_equilibrium(dd, par)
    FUSE.init_core_profiles(dd, par)
    FUSE.init_core_sources(dd, par)

    FUSE.transport_workflow(dd, par, transport_model = :tglfnn, verbose = false)
    run_equilibrium_for_database(dd, par)

    if show_dd_plots
        display(plot(dd.equilibrium, psi_levels_out = [], label = par.general.casename))
        display(plot(dd.core_profiles))
        display(plot(dd.core_sources))
    end

    if !isempty(save_directory)
        IMAS.imas2json(dd, joinpath(save_directory, "$(par.general.casename).json"))
    end
    return dd
end

function validation_workflow(n_samples_per_tokamak::Union{Real,Symbol} = 10; save_directory::String = "", show_dd_plots = false, plot_database = true)

    # Set up the database to run
    run_df = load_dataframe(joinpath(dirname(abspath(@__FILE__)), "..", "..", "sample", "HDB5_compressed.csv"))
    signal_names = ["TOK", "SHOT", "AMIN", "KAPPA", "DELTA", "NEL", "ZEFF", "TAUTH", "RGEO", "BT", "IP", "PNBI", "ENBI", "PICRH", "PECRH", "POHM", "MEFF", "VOL", "AREA", "WTH"]
    run_df = run_df[:, signal_names]
    run_df = run_df[DataFrames.completecases(run_df), :]
    run_df = run_df[(run_df.TOK.!="T10").&(run_df.TOK.!="TDEV").&(run_df.KAPPA.>1.0).&(1.6 .< run_df.MEFF .< 2.2).&(1.1 .< run_df.ZEFF .< 5.9), :]
    if n_samples_per_tokamak !== :all
        tok_list = unique(run_df[:, "TOK"])
        run_df = DataFrames.reduce(
            vcat, [run_df[run_df.TOK.==tok, :][Random.shuffle(1:DataFrames.nrow(run_df[run_df.TOK.==tok, :]))[1:minimum([n_samples_per_tokamak, length(run_df[run_df.TOK.==tok, :][:, "TOK"])])], :] for tok in tok_list]
        )
    end
    # Example of a workflow
    tau_FUSE = Vector{Union{Real,Missing}}(missing, length(run_df[:, "TOK"]))
    tbl = DataFrames.Tables.rowtable(run_df)
    failed_runs_ids = []
    for (idx, data_row) in enumerate(DataFrames.Tables.rows(tbl))
        try
            dd = IMAS.dd()
            par = Parameters(:HDB5; data_row)

            FUSE.init_equilibrium(dd, par)
            FUSE.init_core_profiles(dd, par)
            FUSE.init_core_sources(dd, par)

            FUSE.transport_workflow(dd, par, transport_model = :tglfnn, verbose = false)
            run_equilibrium_for_database(dd, par)

            if show_dd_plots
                display(plot(dd.equilibrium, psi_levels_out = [], label = par.general.casename))
                display(plot(dd.core_profiles))
                display(plot(dd.core_sources))
            end

            if !isempty(save_directory)
                IMAS.imas2json(dd, joinpath(save_directory, "$(par.general.casename).json"))
            end
            tau_FUSE[idx] = @ddtime(dd.summary.global_quantities.tau_energy.value)

        catch e
            push!(failed_runs_ids, idx)
            @show "error = ", e
        end
    end
    run_df[!, :τ_fuse] = tau_FUSE
    println("Number of failed runs $(length(failed_runs_ids)), out of $(length(run_df[:,"TOK"]))")
    if plot_database
        plot_x_y_regression(run_df, "TAUTH", "τ_fuse")
    end
    if !isempty(save_directory)
        CSV.write(joinpath(save_directory, "dataframe.csv"), run_df)
    end
    return run_df
end