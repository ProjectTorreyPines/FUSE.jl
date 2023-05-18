"""
    init(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; do_plot::Bool=false)

Initialize `dd` starting from `ini` and `act` parameters

FUSE provides this high-level `init` function to populate `dd` starting from the `ini` parameters.
This function essentially calls all other `FUSE.init...` functions in FUSE.
For most applications, calling this high level function is sufficient.
"""
function init(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; do_plot::Bool=false)
    TimerOutputs.reset_timer!("init")
    TimerOutputs.@timeit timer "init" begin
        ods_items = []
        # Check what is in the ods to load
        if ini.general.init_from == :ods
            ods_items = keys(IMAS.json2imas(ini.ods.filename))
        end

        # initialize equilibrium
        if !ismissing(ini.equilibrium, :B0) || (:equilibrium ∈ ods_items)
            init_equilibrium(dd, ini, act)
            if do_plot
                display(plot(dd.equilibrium.time_slice[end]))
                plot(dd.equilibrium.time_slice[end]; cx=true, show_x_points=true)
                display(plot!(dd.equilibrium.time_slice[1].boundary, label="Field null"))
            end
        end

        # initialize core profiles
        if !ismissing(ini.core_profiles, :bulk) || (:core_profiles ∈ ods_items)
            init_core_profiles(dd, ini, act)
            if do_plot
                display(plot(dd.core_profiles, legend=:bottomleft))
            end
        end

        # initialize core sources
        if !ismissing(ini.ec_launchers, :power_launched) || !ismissing(ini.ic_antennas, :power_launched) || !ismissing(ini.lh_antennas, :power_launched) || !ismissing(ini.nbi, :power_launched) || (:core_sources ∈ ods_items)
            init_core_sources(dd, ini, act)
            if do_plot
                display(plot(dd.core_sources, legend=:topright))
                display(plot(dd.core_sources, legend=:bottomright; integrated=true))
            end
        end

        # initialize currents
        init_currents(dd, ini, act)

        # initialize build
        if !ismissing(ini.build, :vessel) || !ismissing(ini.build, :layers) || (:build ∈ ods_items)
            init_build(dd, ini, act)
            if do_plot
                plot(dd.equilibrium; cx=true, color=:gray)
                plot!(dd.build)
                display(plot!(dd.build; cx=false))
                display(dd.build.layer)
            end
        end

        # initialize oh and pf coils
        if !ismissing(ini.oh, :n_coils) || (:pf_active ∈ ods_items)
            init_pf_active(dd, ini, act)
            if do_plot
                plot(dd.equilibrium; cx=true, color=:gray)
                plot!(dd.build)
                plot!(dd.build.pf_active.rail)
                display(plot!(dd.pf_active))
            end
        end

        # initialize missing IDSs from ODS (if loading from ODS)
        init_missing_from_ods(dd, ini, act)

        return dd
    end
end

function init(ini::ParametersAllInits, act::ParametersAllActors; do_plot=false)
    dd = IMAS.dd()
    return init(dd, ini, act; do_plot)
end

"""
    init(case::Symbol; do_plot::Bool=false, kw...)

Convenience function to initialize `dd`, `ini`, `act` based on a given use-case.
Returns a tuple with `dd`, `ini`, `act`.
This function is handy if no customization of `ini` or `act` is needed (eg. for regression testing),
otherwise it is recommended to do this in steps:

```julia
ini, act = case_parameters(case::Symbol; kw...)
dd = IMAS.dd()
init(dd, ini, act; do_plot::Bool)
```
"""
function init(case::Symbol; do_plot::Bool=false, kw...)
    ini, act = case_parameters(case; kw...)
    dd = IMAS.dd()
    init(dd, ini, act; do_plot=do_plot)
    return dd, ini, act
end

"""
    consistent_ini_act!(ini::ParametersAllInits, act::ParametersAllActors)

Checks and makes `ini` and `act` consistent with one another
"""
function consistent_ini_act!(ini::ParametersAllInits, act::ParametersAllActors)
    if !ismissing(ini.core_profiles, :T_shaping)
        act.ActorTauenn.T_shaping = ini.core_profiles.T_shaping
    end
    if !ismissing(ini.core_profiles, :T_ratio)
        act.ActorTauenn.T_ratio_pedestal = ini.core_profiles.T_ratio
    end
end