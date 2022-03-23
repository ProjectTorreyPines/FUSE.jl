"""
    init(dd::IMAS.dd, par::Parameters; do_plot = false)

Initialize all IDSs
"""
function init(dd::IMAS.dd, par::Parameters; do_plot = false)
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