function simple_workflow(par::Parameters; do_plot = false)

    dd = IMAS.dd()

    if par.general.init_from == :scalars
        # init equilibrium
        init_equilibrium(dd, par)

        # equilibrium
        eqactor = SolovevEquilibriumActor(dd, symmetric = par.equilibrium.symmetric)
        step(eqactor, verbose = false)
        finalize(eqactor, par.equilibrium.ngrid)

    elseif par.general.init_from == :ods
        dd = IMAS.json2imas(par.ods.filename)
        dd.global_time = dd.equilibrium.time[end]
        IMAS.flux_surfaces(dd.equilibrium)

    elseif par.general.init_from == :gasc
        init_from_gasc(dd, par.gasc.filename, par.gasc.case, par.gasc.no_small_gaps; verbose = true)
    end

    # field null surface
    pushfirst!(dd.equilibrium.time_slice, field_null_surface(dd.equilibrium.time_slice[]))
    pushfirst!(dd.equilibrium.vacuum_toroidal_field.b0, @ddtime(dd.equilibrium.vacuum_toroidal_field.b0))
    pushfirst!(dd.equilibrium.time, -1.0)
    dd.equilibrium.time_slice[1].time = -1.0

    if do_plot
        plot(dd.equilibrium.time_slice[end])
        plot!(dd.equilibrium.time_slice[1].boundary.outline.r, dd.equilibrium.time_slice[1].boundary.outline.z)
    end

    # init radial build
    init_build(dd;
        tf_shape_index = 3,
        is_nuclear_facility = par.build.is_nuclear_facility,
        pf_inside_tf = (par.build.n_pf_coils_inside > 0),
        pf_outside_tf = (par.build.n_pf_coils_outside > 0))

    if do_plot
        plot(dd.equilibrium, color = :gray)
        plot!(dd.build, outline = true)
        display(plot!(dd.build, cx = false))
    end

    # poloidal field coils
    n_coils = [par.build.n_oh_coils]
    if par.build.n_pf_coils_inside > 0
        push!(n_coils, par.build.n_pf_coils_inside)
    end
    if par.build.n_pf_coils_outside > 0
        push!(n_coils, par.build.n_pf_coils_outside)
    end
    pfoptactor = PFcoilsOptActor(dd, n_coils; green_model = par.coil.green_model)
    step(pfoptactor, λ_ψ = 1E-2, λ_null = 1E10, λ_currents = 5E5, λ_strike = 0.0, verbose = true, symmetric = false, maxiter = 1000, optimization_scheme = :rail)
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