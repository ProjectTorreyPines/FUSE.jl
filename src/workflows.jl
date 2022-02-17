function simple_workflow(par::Parameters; do_plot = false)

    dd = IMAS.dd()

    # initialize equilibrium
    init_equilibrium(dd, par)
    if do_plot
        plot(dd.equilibrium.time_slice[end])
        plot!(dd.equilibrium.time_slice[1].boundary.outline.r, dd.equilibrium.time_slice[1].boundary.outline.z)
    end

    # init radial build
    init_build(dd, par)
    if do_plot
        plot(dd.equilibrium, color = :gray)
        plot!(dd.build, outline = true)
        display(plot!(dd.build, cx = false))
    end

    # init oh and pf coils
    init_pf_active(dd, par)

    # optimize coils location
    pfoptactor = PFcoilsOptActor(dd; green_model = par.pf_active.green_model)
    step(pfoptactor, λ_ψ = 1E-2, λ_null = 1E10, λ_currents = 5E5, λ_strike = 0.0, verbose = false, symmetric = false, maxiter = 1000, optimization_scheme = :rail)
    finalize(pfoptactor)

    if do_plot
        display(plot(pfoptactor.trace, :cost))
        display(plot(pfoptactor.trace, :params))

        display(plot(pfoptactor.pf_active, :currents, time_index = 1))
        display(plot(pfoptactor, equilibrium = true, rail = true, time_index = 1))

        display(plot(pfoptactor.pf_active, :currents, time_index = length(dd.equilibrium.time)))
        display(plot(pfoptactor, equilibrium = true, time_index = length(dd.equilibrium.time)))
    end
end