"""
    init(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; do_plot::Bool=false)

Initialize `dd` starting from `ini` and `act` parameters

FUSE provides this high-level `init` function to populate `dd` starting from the `ini` parameters.
This function essentially calls all other `FUSE.init...` functions in FUSE.
For most applications, calling this high level function is sufficient.
"""
function init(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; do_plot::Bool=false, initialize_hardware::Bool=true)
    TimerOutputs.reset_timer!("init")
    TimerOutputs.@timeit timer "init" begin

        # always empty non-hardware IDSs
        empty!(dd.equilibrium)
        empty!(dd.core_profiles)
        empty!(dd.pulse_schedule)
        empty!(dd.core_sources)
        empty!(dd.summary)

        # set the dd.global time to when simulation starts
        dd.global_time = ini.time.simulation_start

        # we make a copy because we overwrite some of the parameters internally
        # and we want this function to work always the same for subsequent calls
        ini = deepcopy(ini)
        act = deepcopy(act)

        # load ods once if needed
        dd1 = ini_from_ods!(ini)

        # Makes `ini` and `act` self-consistent and consistent with one another
        consistent_ini_act!(ini, act)

        # initialize pulse_schedule
        if !initialize_hardware || !ismissing(ini.equilibrium, :B0) || !isempty(dd1.equilibrium) || !isempty(dd1.pulse_schedule)
            init_pulse_schedule!(dd, ini, act, dd1)
            if do_plot
                display(plot(dd.pulse_schedule))
            end
        end

        # initialize equilibrium
        if !initialize_hardware || !ismissing(ini.equilibrium, :B0) || !isempty(dd1.equilibrium)
            if ini.general.init_from == :ods && !isempty(dd1.pf_active.coil)
                init_pf_active!(dd, ini, act, dd1)
            end
            for coil in dd.pf_active.coil
                empty!(coil.current)
            end
            init_equilibrium!(dd, ini, act, dd1)
            if do_plot
                display(plot(dd.equilibrium.time_slice[end]))
                plot(dd.equilibrium.time_slice[end]; cx=true, show_x_points=true)
                display(plot!(dd.equilibrium.time_slice[1].boundary; label="Field null"))
            end
        end

        # initialize core profiles
        if !initialize_hardware || !ismissing(ini.core_profiles, :bulk) || !isempty(dd1.core_profiles)
            init_core_profiles!(dd, ini, act, dd1)
            if do_plot
                display(plot(dd.core_profiles; legend=:bottomleft))
            end
        end

        # initialize core sources
        if !initialize_hardware || !isempty(ini.ec_launcher) || !isempty(ini.ic_antenna) || !isempty(ini.lh_antenna) || !isempty(ini.nb_unit) || !isempty(dd1.core_sources)
            init_core_sources!(dd, ini, act, dd1)
            if do_plot
                display(plot(dd.core_sources; legend=:topright))
                display(plot(dd.core_sources; legend=:bottomright, integrated=true))
            end
        end

        # initialize currents
        init_currents!(dd, ini, act, dd1)

        # initialize build
        if initialize_hardware && (!isempty(ini.build.layers) || !isempty(dd1.build))
            init_build!(dd, ini, act, dd1)
            if do_plot
                plot(dd.equilibrium; cx=true, color=:gray)
                plot!(dd.build)
                display(plot!(dd.build; cx=false))
                display(dd.build.layer)
            end
        end

        # initialize oh and pf coils
        if initialize_hardware && (!ismissing(ini.oh, :n_coils) || !isempty(dd1.pf_active.coil))
            init_pf_active!(dd, ini, act, dd1)
            if do_plot
                plot(dd.equilibrium; cx=true, color=:gray)
                plot!(dd.build)
                plot!(dd.build.pf_active.rail)
                display(plot!(dd.pf_active))
            end
        end

        # initialize requirements
        init_requirements!(dd, ini, act, dd1)

        # initialize missing IDSs from ODS (if loading from ODS)
        init_missing_from_ods!(dd, ini, act, dd1)

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

Makes `ini` and `act` self-consistent and consistent with one another

NOTE: operates in place
"""
function consistent_ini_act!(ini::ParametersAllInits, act::ParametersAllActors)

    if !ismissing(ini.core_profiles, :T_ratio)
        act.ActorFixedProfiles.T_ratio_core = ini.core_profiles.T_ratio
        act.ActorPedestal.T_ratio_pedestal = ini.core_profiles.T_ratio
        act.ActorTauenn.T_ratio_pedestal = ini.core_profiles.T_ratio
    end

    if !ismissing(ini.core_profiles, :T_shaping)
        act.ActorFixedProfiles.T_shaping = ini.core_profiles.T_shaping
        act.ActorTauenn.T_shaping = ini.core_profiles.T_shaping
    end

    if !ismissing(ini.core_profiles, :n_shaping)
        act.ActorFixedProfiles.n_shaping = ini.core_profiles.n_shaping
    end

    if !ismissing(ini.equilibrium, :xpoints)
        act.ActorPFdesign.symmetric = ini.equilibrium.xpoints in [:none, :double]
    end
end
