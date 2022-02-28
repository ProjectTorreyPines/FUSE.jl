function init_workflow(dd::IMAS.dd, par::Parameters; do_plot = false)
    # initialize equilibrium
    init_equilibrium(dd, par)
    if do_plot
        plot(dd.equilibrium.time_slice[end])
        plot!(dd.equilibrium.time_slice[1].boundary.outline.r, dd.equilibrium.time_slice[1].boundary.outline.z)
    end

    # initialize build
    init_build(dd, par)
    if do_plot
        plot(dd.equilibrium, color = :gray)
        plot!(dd.build)
        display(plot!(dd.build, cx = false))
    end

    # initialize oh and pf coils
    init_pf_active(dd, par)
    if do_plot
        plot(dd.equilibrium, color = :gray)
        plot!(dd.build)
        plot!(dd.build.pf_active.rail)
        display(plot!(dd.pf_active))
    end

    # initialize core profiles
    init_core_profiles(dd, par)
    if do_plot
        display(plot(dd.core_profiles))
    end

    # initialize core sources
    init_core_sources(dd, par)
    if do_plot
        display(plot(dd.core_sources))
        display(plot(dd.core_sources; integrated = true))
    end

    # initialize missing IDSs (if loading from ODS)
    init_missing(dd, par)

    return dd
end

function pf_optim_workflow(dd::IMAS.dd, par::Parameters; update_eq_in = false, do_plot = false)
    # initialize actor
    actor = PFcoilsOptActor(dd; green_model = par.pf_active.green_model, symmetric = par.build.symmetric)

    # optimize coil location only considering equilibria (disregard field-null)
    step(actor, λ_ψ = 1E-2, λ_null = 1E10, λ_currents = 1E6, λ_strike = 0.0, verbose = false, maxiter = 1000, optimization_scheme = :rail)
    finalize(actor)

    if do_plot
        display(plot(actor.trace, :cost))
        display(plot(actor.trace, :params))
    end

    # find coil currents for both field-null and equilibria
    step(actor, λ_ψ = 1E-2, λ_null = 1E-2, λ_currents = 1E6, λ_strike = 0.0, verbose = false, maxiter = 1000, optimization_scheme = :static)
    finalize(actor; update_eq_in)

    if do_plot
        display(plot(actor.pf_active, :currents, time_index = 1))
        display(plot(actor, equilibrium = true, rail = true, time_index = 1))

        display(plot(actor.pf_active, :currents, time_index = length(dd.equilibrium.time)))
        display(plot(actor, equilibrium = true, time_index = length(dd.equilibrium.time)))
    end

    return dd
end

function build_workflow(dd::IMAS.dd, par::Parameters; rebuild_wall=(par.general.init_from != :ods), do_plot = false)
    # optimize pf_active and update equilibrium
    pf_optim_workflow(dd::IMAS.dd, par::Parameters; update_eq_in = true, do_plot)
    
    if rebuild_wall
        # regenerate build based on new equilibrium
        empty!(dd.wall)
        init_build(dd, par)

        # re-optimize pf_active
        pf_optim_workflow(dd::IMAS.dd, par::Parameters; update_eq_in = false, do_plot=false)
    end

    if do_plot
        plot(dd.equilibrium, color = :gray)
        plot!(dd.build)
        plot!(dd.build.pf_active.rail)
        display(plot!(dd.pf_active))
    end
end

function transport_workflow(dd::IMAS.dd, par::Parameters; do_plot = false, verbose = false, kw...)
    if do_plot
        plot(dd.core_profiles; color = :gray, label = "")
    end
    actor = TaueNNactor(dd; kw...)
    step(actor; verbose = verbose)
    finalize(actor)
    if do_plot
        display(plot!(dd.core_profiles))
    end
    return dd
end