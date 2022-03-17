using DataFrames
using CSV
using Random
function init_dataframe_HDB5()
    # For description of variables see https://osf.io/593q6/
    filename = joinpath(dirname(abspath(@__FILE__)), "..", "..", "sample", "HDB5_compressed.csv")
    return CSV.read(filename, DataFrame)
end

function run_equilibrium_for_database(dd, case)
    # Set beta_normal from equilbrium to the kinetic beta_n
    if !isempty(dd.core_profiles.profiles_1d)
        dd.equilibrium.time_slice[].global_quantities.beta_normal = @ddtime(dd.summary.global_quantities.beta_tor_thermal_norm.value)
    end
    eqactor = FUSE.SolovevEquilibriumActor(dd.equilibrium)
    FUSE.step(eqactor, verbose = false)
    FUSE.finalize(eqactor)

    # correct equilibrium volume and area
    eqt1d = dd.equilibrium.time_slice[].profiles_1d
    eqt1d.volume .*= case[:VOL] / eqt1d.volume[end]
    eqt1d.area .*= case[:AREA] / eqt1d.area[end]
end

"""
    function validation_workflow(
        n_samples::Union{Real,Symbol},
        save_directory::String,
        show_equilibrium::Bool,
        plot_database::Bool)

Runs n_samples of the HDB5 database and stores results in save_directory
"""

function validation_workflow(n_samples::Union{Real,Symbol} = 10; save_directory::String = "", show_dd_plots = false, plot_database = true)

    # Set up the database to run
    run_df = init_dataframe_HDB5()
    signal_names = ["TOK", "SHOT", "AMIN", "KAPPA", "DELTA", "NEL", "ZEFF", "TAUTH", "RGEO", "BT", "IP", "PNBI", "ENBI", "PICRH", "PECRH", "POHM", "MEFF", "VOL", "AREA", "WTH"]
    run_df = run_df[:, signal_names]
    run_df = run_df[completecases(run_df), :]
    run_df = run_df[(run_df.TOK.!="T10").&(run_df.TOK.!="TDEV").&(run_df.KAPPA.>1.0).&(1.6 .< run_df.MEFF .< 2.2).&(1.1 .< run_df.ZEFF .< 5.9), :]
    if n_samples !== :all
        run_df = run_df[shuffle(1:nrow(run_df))[1:n_samples], :]
    end
    # Example of a workflow
    tau_FUSE = Vector{Union{Real,Missing}}(missing, length(run_df[:, "TOK"]))
    tbl = Tables.rowtable(run_df)
    failed_runs_ids = []
    for (idx, case) in enumerate(Tables.rows(tbl))
        try
            dd = IMAS.dd()
            par = Parameters(:HDB5, case = case)
            case_name = "$(case[:TOK])_$(case[:SHOT])"

            FUSE.init_equilibrium(dd, par)

            run_equilibrium_for_database(dd, case)

            FUSE.init_core_profiles(dd, par)
            FUSE.init_core_sources(dd, par)

            FUSE.transport_workflow(dd, par, transport_model = :tglfnn, verbose = false)
            run_equilibrium_for_database(dd, case)

            if show_dd_plots
                display(plot(dd.equilibrium, psi_levels_out = [], label = case_name))
                display(plot(dd.core_profiles))
                display(plot(dd.core_sources))
            end

            if !isempty(save_directory)
                IMAS.imas2json(dd, joinpath(save_directory, "$case_name.json"))
            end
            tau_FUSE[idx] = @ddtime(dd.summary.global_quantities.tau_energy.value)

        catch e
            push!(failed_runs_ids, idx)
            @show "error = ", e
        end
    end
    run_df[!, :τ_fuse] = tau_FUSE
    if plot_database
        x_ylim = [1e-2, 10.]
        println("Number of failed runs $(length(failed_runs_ids)), out of $(length(run_df[:,"TOK"]))")
        plot(run_df[:, :TAUTH], run_df[:, :τ_fuse], seriestype = :scatter, xaxis = :log, yaxis = :log, ylim = x_ylim, xlim = x_ylim, xlabel = "τ_e_exp [s]", ylabel = "τ_e_FUSE [s]", label = nothing)
        plot!([0.5 * x_ylim[1], 0.5 * x_ylim[2]], [2 * x_ylim[1], 2 * x_ylim[2]], linestyle = :dash, label = "+50%")
        plot!([2 * x_ylim[1], 2 * x_ylim[2]], [0.5 * x_ylim[1], 0.5 * x_ylim[2]], linestyle = :dash, label = "-50%", legend = :topleft)
        display(plot!([x_ylim[1], x_ylim[2]], [x_ylim[1], x_ylim[2]], label = nothing))
    end
    if !isempty(save_directory)
        CSV.write(joinpath(save_directory, "dataframe.csv"), run_df)
    end
    return run_df
end