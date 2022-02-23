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
        plot!(dd.build, outline = true)
        display(plot!(dd.build, cx = false))
    end

    # initialize oh and pf coils
    init_pf_active(dd, par)
    if do_plot
        plot(dd.equilibrium, color = :gray)
        plot!(dd.build, outline = true)
        plot!(dd.build.pf_coils_rail)
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
    end

    return dd
end

function pf_optim_workflow(dd::IMAS.dd, par::Parameters; do_plot = false)
    pfoptactor = PFcoilsOptActor(dd; green_model = par.pf_active.green_model)
    step(pfoptactor, λ_ψ = 1E-2, λ_null = 1E10, λ_currents = 1E6, λ_strike = 0.0, verbose = false, symmetric = false, maxiter = 1000, optimization_scheme = :rail)
    finalize(pfoptactor)

    if do_plot
        display(plot(pfoptactor.trace, :cost))
        display(plot(pfoptactor.trace, :params))

    end

    step(pfoptactor, λ_ψ = 1E-2, λ_null = 1E-2, λ_currents = 1E6, λ_strike = 0.0, verbose = false, symmetric = false, maxiter = 1000, optimization_scheme = :static)
    finalize(pfoptactor)

    if do_plot
        display(plot(pfoptactor.pf_active, :currents, time_index = 1))
        display(plot(pfoptactor, equilibrium = true, rail = true, time_index = 1))

        display(plot(pfoptactor.pf_active, :currents, time_index = length(dd.equilibrium.time)))
        display(plot(pfoptactor, equilibrium = true, time_index = length(dd.equilibrium.time)))
    end

    return dd
end

function transport_workflow(dd::IMAS.dd, par::Parameters; do_plot = false)
    if do_plot
        plot(dd.core_profiles; color = :gray, label="")
    end
    tauennactor = TaueNNactor(dd, use_tglfnn = false, error = 1E-2)
    step(tauennactor; verbose = true)
    finalize(tauennactor)
    if do_plot
        display(plot!(dd.core_profiles))
    end
    return dd
end