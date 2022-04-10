"""
    init(dd::IMAS.dd, ini::InitParameters, act::ActorParameters; do_plot=false)

Initialize all IDSs if there are parameters for it or is initialized from ods
"""
function init(dd::IMAS.dd, ini::InitParameters, act::ActorParameters; do_plot=false)
    ods_items = []
    # Check what is in the ods to load
    if ini.general.init_from == :ods
        ods_items = keys(IMAS.json2imas(ini.ods.filename))
    end

    # initialize equilibrium
    if !ismissing(ini.equilibrium, :B0) || :equilibrium ∈ ods_items
        init_equilibrium(dd, ini, act)
        if do_plot
            plot(dd.equilibrium.time_slice[end]; x_point=true)
            display(plot!(dd.equilibrium.time_slice[1].boundary.outline.r, dd.equilibrium.time_slice[1].boundary.outline.z,label="Field null"))
        end
    end

    # initialize build
    if !ismissing(ini.build, :vessel) || !ismissing(ini.build, :layers) || :build ∈ ods_items
        init_build(dd, ini, act)
        if do_plot
            plot(dd.equilibrium, color=:gray)
            plot!(dd.build)
            display(plot!(dd.build, cx=false))
        end
    end

    # initialize oh and pf coils
    if !ismissing(ini.pf_active, :n_oh_coils) || :pf_active ∈ ods_items
        init_pf_active(dd, ini, act)
        if do_plot
            plot(dd.equilibrium, color=:gray)
            plot!(dd.build)
            plot!(dd.build.pf_active.rail)
            display(plot!(dd.pf_active))
        end
    end

    # initialize core profiles
    if !ismissing(ini.core_profiles, :bulk) || :core_profiles ∈ ods_items
        init_core_profiles(dd, ini, act)
        if do_plot
            display(plot(dd.core_profiles))
        end
    end

    # initialize core sources
    if !ismissing(ini.ec, :power_launched) || !ismissing(ini.ic, :power_launched) || !ismissing(ini.lh, :power_launched) || !ismissing(ini.nbi, :power_launched) || :core_sources ∈ ods_items
        init_core_sources(dd, ini, act)
        if do_plot
            display(plot(dd.core_sources))
            display(plot(dd.core_sources; integrated=true))
        end
    end

    # initialize missing IDSs (if loading from ODS)
    init_missing(dd, ini, act)

    return dd
end