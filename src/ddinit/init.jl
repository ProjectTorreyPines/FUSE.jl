"""
    init(dd::IMAS.dd, par::Parameters; do_plot = false)

Initialize all IDSs
"""
function init(dd::IMAS.dd, par::Parameters; do_plot=false)
    # initialize equilibrium
    init_equilibrium(dd, par)
    if do_plot
        plot(dd.equilibrium.time_slice[end])
        display(plot!(dd.equilibrium.time_slice[1].boundary.outline.r, dd.equilibrium.time_slice[1].boundary.outline.z))
    end

    # initialize build
    if !ismissing(par.build, :vessel) || !ismissing(par.build, :layers)
        init_build(dd, par)
        if do_plot
            plot(dd.equilibrium, color=:gray)
            plot!(dd.build)
            display(plot!(dd.build, cx=false))
        end
    end

    # initialize oh and pf coils
    if !ismissing(par.pf_active, :n_oh_coils)
        init_pf_active(dd, par)
        if do_plot
            plot(dd.equilibrium, color=:gray)
            plot!(dd.build)
            plot!(dd.build.pf_active.rail)
            display(plot!(dd.pf_active))
        end
    end

    # initialize core profiles
    if !ismissing(par.core_profiles, :bulk)
        init_core_profiles(dd, par)
        if do_plot
            display(plot(dd.core_profiles))
        end
    end

    # initialize core sources
    if !ismissing(par.ec, :power_launched) || !ismissing(par.ic, :power_launched) || !ismissing(par.lh, :power_launched) || !ismissing(par.nbi, :power_launched)
        init_core_sources(dd, par)
        if do_plot
            display(plot(dd.core_sources))
            display(plot(dd.core_sources; integrated=true))
        end
    end

    # initialize missing IDSs (if loading from ODS)
    init_missing(dd, par)

    return dd
end