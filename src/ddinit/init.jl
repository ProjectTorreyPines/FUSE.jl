"""
    init(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; kw...)

Like `init!(dd, ini, act)` function, but does not modify `ini` and `act`
"""
function init(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; kw...)
    ini = deepcopy(ini)
    act = deepcopy(act)
    return init!(dd, ini, act; kw...)
end

"""
    init!(
        dd::IMAS.dd,
        ini::ParametersAllInits,
        act::ParametersAllActors;
        do_plot::Bool=false,
        initialize_hardware::Bool=true,
        initialize_pulse_schedule::Bool=true,
        verbose::Bool=false)

Initialize `dd` starting from `ini` and `act` parameters

FUSE provides this high-level `init` function to populate `dd` starting from the `ini` parameters

This function calls all other `FUSE.init...` functions in FUSE

Modifies `ini` and `act` in place
"""
function init!(
    dd::IMAS.dd{T},
    ini::ParametersAllInits,
    act::ParametersAllActors;
    do_plot::Bool=false,
    initialize_hardware::Bool=true,
    initialize_pulse_schedule::Bool=true,
    restore_expressions::Bool=true,
    verbose::Bool=false) where {T<:Real}

    TimerOutputs.reset_timer!("init")
    TimerOutputs.@timeit timer "init" begin

        # always empty non-hardware IDSs
        empty!(dd.equilibrium)
        empty!(dd.core_profiles)
        empty!(dd.edge_profiles)
        empty!(dd.core_sources)
        empty!(dd.summary)

        # optionally re-initialize pulse_schedule
        if initialize_pulse_schedule
            empty!(dd.pulse_schedule)
        end

        # set the dd.global time to when simulation starts
        dd.global_time = ini.time.simulation_start

        # load ods once
        verbose && @info "INIT: load_ods"
        dd1 = set_ini_act_from_ods!(ini, act)

        # Here we delete fields from the ODS for which we know FUSE has expressions for.
        # Besides ensuring consistency, this is done because some FUSE workflows in fact expect certain fields to be expressions!
        if restore_expressions
            FUSE.restore_init_expressions!(dd1; verbose=false)
        end

        # initialize pulse_schedule
        ps_was_set = false
        if ismissing(dd.pulse_schedule.flux_control, :time) || isempty(dd.pulse_schedule.flux_control.time)
            empty!(dd.pulse_schedule)
            verbose && @info "INIT: init_pulse_schedule"
            init_pulse_schedule!(dd, ini, act, dd1)
            if do_plot
                display(plot(dd.pulse_schedule))
            end
            ps_was_set = true
        end

        # wall
        if ini.general.init_from == :ods
            if !isempty(dd1.wall.description_2d)
                dd.wall = deepcopy(dd1.wall)
            end
        end

        # initialize equilibrium
        if !initialize_hardware || !ismissing(ini.equilibrium, :B0) || !isempty(dd1.equilibrium)
            if ini.general.init_from == :ods && !isempty(dd1.pf_active.coil)
                verbose && @info "INIT: init_pf_active"
                init_pf_active!(dd, ini, act, dd1)
            end
            if isempty(dd.pf_active.coil)
                # add some coils encircling the plasma so that the equilibrium solvers can get free-boundary solution
                bnd_r, bnd_z = IMAS.boundary(dd.pulse_schedule.position_control)
                dd.pf_active.coil = encircling_coils(bnd_r, bnd_z, sum(extrema(bnd_r)) / 2.0, sum(extrema(bnd_z)) / 2.0, 8)
                act.ActorCXbuild.layers_aware_of_pf_coils = false
            else
                for coil in dd.pf_active.coil
                    empty!(coil.current)
                end
            end
            verbose && @info "INIT: init_equilibrium"
            init_equilibrium!(dd, ini, act, dd1)
            if do_plot
                display(plot(dd.equilibrium.time_slice[end]))
                plot(dd.equilibrium.time_slice[end]; cx=true, show_x_points=true)
                display(plot!(dd.equilibrium.time_slice[1].boundary; label="Field null"))
            end
        end

        # initialize core profiles
        if !initialize_hardware || !ismissing(ini.core_profiles, :bulk) || !isempty(dd1.core_profiles)
            verbose && @info "INIT: init_core_profiles"
            init_core_profiles!(dd, ini, act, dd1)
            if do_plot
                display(plot(dd.core_profiles; legend=:bottomleft))
            end
        end

        # initialize edge profiles
        if !initialize_hardware || !ismissing(ini.core_profiles, :bulk) || !isempty(dd1.edge_profiles)
            verbose && @info "INIT: init_edge_profiles"
            init_edge_profiles!(dd, ini, act, dd1)
            if do_plot
                display(plot(dd.edge_profiles; legend=:bottomleft))
            end
        end

        # initialize build
        if initialize_hardware && (!isempty(ini.build.layers) || !isempty(dd1.build))
            verbose && @info "INIT: init_build"
            init_build!(dd, ini, act, dd1)
            if do_plot
                plot(dd.equilibrium; cx=true, color=:gray)
                plot!(dd.build; equilibrium=false, pf_active=false)
                display(plot!(dd.build; cx=false))
                display(dd.build.layer)
            end
        end

        # initialize core sources and HCD
        if (
            !initialize_hardware || !isempty(ini.ec_launcher) || !isempty(ini.pellet_launcher) || !isempty(ini.ic_antenna) || !isempty(ini.lh_antenna) || !isempty(ini.nb_unit) ||
            !isempty(dd1.core_sources)
        )
            verbose && @info "INIT: init_core_sources and HCD"
            if ismissing(ini.hcd, :power_scaling_cost_function)
                init_core_sources!(dd, ini, act, dd1)
                init_hcd!(dd, ini, act, dd1)
            else
                # get an estimate for Ohmic heating before scaling HCD powers
                ini0 = deepcopy(ini)
                ps = dd.pulse_schedule
                ps0 = deepcopy(ps)
                function scale_power_tau_cost(scale; dd, ps, ini, ps0, ini0, power_scaling_cost_function)
                    scale_powers!(ps, ini, ps0, ini0, scale)
                    init_core_sources!(dd, ini, act, dd1)
                    init_hcd!(dd, ini, act, dd1)
                    init_currents!(dd, ini, act, dd1) # to pick up ohmic power
                    error = abs(power_scaling_cost_function(dd))
                    return error
                end

                try
                    old_logging = actor_logging(dd, false)
                    res =
                        Optim.optimize(
                            scale -> scale_power_tau_cost(scale; dd, ps, ini, ps0, ini0, ini.hcd.power_scaling_cost_function),
                            0.0,
                            100,
                            Optim.Brent();
                            abs_tol=1E-3
                        )
                    actor_logging(dd, old_logging)
                    scale_power_tau_cost(res.minimizer; dd, ps, ini, ps0, ini0, ini.hcd.power_scaling_cost_function)
                catch e
                    actor_logging(dd, old_logging)
                    rethrow(e)
                end

            end

            if do_plot
                display(plot(dd.core_sources; legend=:topright))
                display(plot(dd.core_sources; legend=:bottomright, integrated=true))
            end
        end

        # initialize currents
        verbose && @info "INIT: init_currents"
        init_currents!(dd, ini, act, dd1)

        # initialize oh and pf coils
        n_coils = [length(layer.coils_inside) for layer in dd.build.layer if !ismissing(layer, :coils_inside)]
        if initialize_hardware && !isempty(n_coils)
            verbose && @info "INIT: init_pf_active"
            init_pf_active!(dd, ini, act, dd1)
            if do_plot
                plot(dd.equilibrium; cx=true, color=:gray)
                plot!(dd.build; equilibrium=false, pf_active=false)
                plot!(dd.build.pf_active.rail)
                display(plot!(dd.pf_active))
            end
        end

        # pf_passive
        if !isempty(dd.build.layer)
            if ini.general.init_from == :ods && !isempty(dd1.pf_passive.loop)
                dd.pf_passive = deepcopy(dd1.pf_passive)
            else
                FUSE.ActorPassiveStructures(dd, act)
            end
        end

        # initialize balance of plant
        verbose && @info "INIT: init_balance_of_plant"
        init_balance_of_plant!(dd, ini, act, dd1)

        # initialize requirements
        verbose && @info "INIT: init_requirements"
        init_requirements!(dd, ini, act, dd1)

        # initialize missing IDSs from ODS (if loading from ODS)
        verbose && @info "INIT: init_missing_from_ods"
        init_missing_from_ods!(dd, ini, act, dd1)

        # add strike point information to pulse_schedule
        fw = IMAS.first_wall(dd.wall)
        if !isempty(fw.r) && ps_was_set
            RXX = Vector{T}[]
            ZXX = Vector{T}[]
            for eqt in dd.equilibrium.time_slice
                if eqt.global_quantities.ip == 0.0
                    Rxx, Zxx = T[], T[]
                else
                    psi_boundaries = (last_closed=eqt.boundary.psi, first_open=eqt.boundary_separatrix.psi)
                    Rxx, Zxx, _ = IMAS.find_strike_points(eqt, fw.r, fw.z, psi_boundaries.last_closed, psi_boundaries.first_open)
                end
                push!(RXX, Rxx)
                push!(ZXX, Zxx)
            end
            N = maximum(collect(map(length, RXX)))
            pc = dd.pulse_schedule.position_control
            resize!(pc.strike_point, N)
            for k in 1:N
                pc.strike_point[k].r.reference = IMAS.interp1d(dd.equilibrium.time, [k <= length(Rxx) ? Rxx[k] : T(NaN) for Rxx in RXX], :constant).(pc.time)
                pc.strike_point[k].z.reference = IMAS.interp1d(dd.equilibrium.time, [k <= length(Zxx) ? Zxx[k] : T(NaN) for Zxx in ZXX], :constant).(pc.time)
            end
        end

        # trim core_profiles and edge_profiles data before the first equilibrium since things are really not robust against that
        # also trim other IDSs not to go past equilibrium.time[end]
        if dd.equilibrium.time[1] != dd.equilibrium.time[end]
            IMAS.trim_time!(dd, (-Inf, dd.equilibrium.time[end]))
            IMAS.trim_time!(dd.core_profiles, (dd.equilibrium.time[1], dd.equilibrium.time[end]))
            IMAS.trim_time!(dd.edge_profiles, (dd.equilibrium.time[1], dd.equilibrium.time[end]))
        end

        # setup ActorReplay
        act.ActorReplay.replay_dd = dd1

        return dd
    end
end

function scale_powers!(ps::IMAS.pulse_schedule, ini::ParametersAllInits, ps0::IMAS.pulse_schedule, ini0::ParametersAllInits, scale::Float64)
    # pulse_shedule
    for (beam, beam0) in zip(ps.ec.beam, ps0.ec.beam)
        beam.power_launched.reference .= beam0.power_launched.reference .* scale
    end
    for (antenna, antenna0) in zip(ps.ic.antenna, ps0.ic.antenna)
        antenna.power.reference .= antenna0.power.reference .* scale
    end
    for (antenna, antenna0) in zip(ps.lh.antenna, ps0.lh.antenna)
        antenna.power.reference .= antenna0.power.reference .* scale
    end
    for (unit, unit0) in zip(ps.nbi.unit, ps0.nbi.unit)
        unit.power.reference .= unit0.power.reference .* scale
    end

    #ini
    for (launcher, launcher0) in zip(ini.ec_launcher, ini0.ec_launcher)
        launcher.power_launched = launcher0.power_launched .* scale
    end
    for (antenna, antenna0) in zip(ini.ic_antenna, ini0.ic_antenna)
        antenna.power_launched = antenna0.power_launched .* scale
    end
    for (antenna, antenna0) in zip(ini.lh_antenna, ini0.lh_antenna)
        antenna.power_launched = antenna0.power_launched .* scale
    end
    for (unit, unit0) in zip(ini.nb_unit, ini0.nb_unit)
        unit.power_launched = unit0.power_launched .* scale
    end
end

function init(ini::ParametersAllInits, act::ParametersAllActors; do_plot=false)
    dd = IMAS.dd()
    return init(dd, ini, act; do_plot)
end

"""
    init(case::Symbol; do_plot::Bool=false, kw...)

Initialize `dd`, `ini`, `act` based on a given use-case.

Returns a tuple with `dd`, `ini`, `act`.
"""
function init(case::Symbol; do_plot::Bool=false, kw...)
    ini, act = case_parameters(case; kw...)
    dd = IMAS.dd()
    init(dd, ini, act; do_plot=do_plot)
    return dd, ini, act
end
