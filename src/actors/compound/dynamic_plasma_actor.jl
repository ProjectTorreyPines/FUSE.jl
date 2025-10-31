#= ================== =#
#  ActorDynamicPlasma  #
#= ================== =#
@actor_parameters_struct ActorDynamicPlasma{T} begin
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution")
    evolve_transport::Entry{Bool} = Entry{Bool}("-", "Evolve the transport"; default=true)
    evolve_pedestal::Entry{Bool} = Entry{Bool}("-", "Evolve the pedestal"; default=true)
    evolve_hcd::Entry{Bool} = Entry{Bool}("-", "Evolve the heating and current drive"; default=true)
    evolve_current::Entry{Bool} = Entry{Bool}("-", "Evolve the plasma current"; default=true)
    evolve_equilibrium::Entry{Bool} = Entry{Bool}("-", "Evolve the equilibrium"; default=true)
    evolve_pf_active::Entry{Bool} = Entry{Bool}("-", "Evolve the PF currents"; default=false)
    ip_controller::Entry{Bool} = Entry{Bool}("-", "Use controller to change v_loop to match desired Ip"; default=false)
    time_derivatives_sources::Entry{Bool} = Entry{Bool}("-", "Include time-derivative sources"; default=true)
    #== display and debugging parameters ==#
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorDynamicPlasma{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorDynamicPlasma{P}}
    act::ParametersAllActors{P}
    actor_tr::ActorCoreTransport{D,P}
    actor_ped::ActorPedestal{D,P}
    actor_hc::ActorHCD{D,P}
    actor_jt::ActorCurrent{D,P}
    actor_eq::ActorEquilibrium{D,P}
    actor_pf::ActorPFactive{D,P}
    actor_saw::ActorSawteeth{D,P}
end

"""
    ActorDynamicPlasma(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves plasma profiles forward in time using coupled physics models.

This compound actor advances the plasma state through time by coordinating multiple
physics models with proper temporal splitting and time step management. It handles
the complex coupling between transport, current evolution, heating, equilibrium,
and magnetic control systems.

Time integration strategy:

  - Uses second-order accurate time-stepping with sub-cycling
  - Splits fast (equilibrium, PF) and slow (transport, pedestal) physics appropriately
  - Handles time-dependent boundary conditions from pulse schedule
  - Supports plasma current control via feedback systems

Evolution workflow (each time step):

 1. **Phase 1**: Current evolution and heating sources (δt steps)
 2. **Phase 2**: Pedestal and transport evolution (δt steps)
 3. **Both phases**: Sawteeth, equilibrium, and PF control (δt/2 sub-steps)

Key features:

  - Configurable evolution of individual physics components
  - Optional plasma current feedback control via loop voltage
  - Time derivative sources for energy and particle balance
  - Progress tracking and convergence monitoring
  - Automatic cleanup of post-simulation time data

Physics models coordinated:

  - **ActorCurrent**: Current density diffusion and resistive evolution
  - **ActorHCD**: Time-dependent heating and current drive
  - **ActorPedestal**: Pedestal evolution and L-H transitions
  - **ActorCoreTransport**: Transport flux evolution
  - **ActorSawteeth**: Sawtooth instability cycling
  - **ActorEquilibrium**: MHD equilibrium with evolving profiles
  - **ActorPFactive**: PF coil current control for plasma shape

Control options:

  - `ip_controller`: Feedback control of plasma current via loop voltage
  - Individual physics component enable/disable switches
  - Configurable time step size and number of steps
"""
function ActorDynamicPlasma(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorDynamicPlasma(dd, act.ActorDynamicPlasma, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorDynamicPlasma(dd::IMAS.dd, par::FUSEparameters__ActorDynamicPlasma, act::ParametersAllActors; kw...)
    logging_actor_init(ActorDynamicPlasma)
    par = OverrideParameters(par; kw...)

    actor_tr = ActorCoreTransport(dd, act.ActorCoreTransport, act)

    # allows users to hardwire `rho_nml` and `rho_ped` (same logic here as in ActorStationaryPlasma)
    if act.ActorCoreTransport.model == :FluxMatcher && ismissing(act.ActorPedestal, :rho_nml)
        rho_nml = actor_tr.tr_actor.par.rho_transport[end-1]
    elseif ismissing(act.ActorPedestal, :rho_nml)
        rho_nml = act.ActorFluxMatcher.rho_transport[end-1]
    else
        rho_nml = act.ActorPedestal.rho_nml
    end
    if act.ActorCoreTransport.model == :FluxMatcher && ismissing(act.ActorPedestal, :rho_ped)
        rho_ped = actor_tr.tr_actor.par.rho_transport[end]
    elseif ismissing(act.ActorPedestal, :rho_ped)
        rho_ped = act.ActorFluxMatcher.rho_transport[end]
    else
        rho_ped = act.ActorPedestal.rho_ped
    end
    actor_ped = ActorPedestal(
        dd,
        act.ActorPedestal,
        act;
        ip_from=:pulse_schedule,
        βn_from=:core_profiles,
        ne_from=:pulse_schedule,
        zeff_from=:pulse_schedule,
        rho_nml,
        rho_ped)

    actor_hc = ActorHCD(dd, act.ActorHCD, act)

    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; ip_from=:pulse_schedule, vloop_from=:pulse_schedule)

    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:pulse_schedule)

    actor_pf = ActorPFactive(dd, act.ActorPFactive)

    actor_saw = ActorSawteeth(dd, act.ActorSawteeth)

    return ActorDynamicPlasma(dd, par, act, actor_tr, actor_ped, actor_hc, actor_jt, actor_eq, actor_pf, actor_saw)
end

function _step(actor::ActorDynamicPlasma)
    dd = actor.dd
    par = actor.par

    Δt = par.Δt
    Nt = par.Nt
    δt = Δt / Nt

    t0 = dd.global_time
    t1 = dd.global_time + Δt

    # remove time dependent data after global_time
    IMAS.trim_time!(actor.dd, (-Inf, dd.global_time); trim_pulse_schedule=false)

    substeps_per_2loop = 9
    ProgressMeter.ijulia_behavior(:clear)
    prog = ProgressMeter.Progress(Nt * substeps_per_2loop; dt=0.0, showspeed=true, enabled=par.verbose)
    old_logging = actor_logging(dd, false)

    try
        # NOTE: δt is a full step, some actors are called every 1/2 step
        for (kk, time0) in enumerate(range(t0, t1, 2 * Nt + 1)[2:end])
            phase = mod(kk + 1, 2) + 1 # phase can be either 1 or 2, we start with 1
            progr = (prog, t0, t1, phase)

            # Prepare time dependent arrays of structures
            # NOTE: dd.core_profiles is different because it is updated
            # by actor_jt at the 1/2 steps, but also
            # by actor_tr and actor_ped at the 1/2 steps.
            # For dd.core_profiles we thus create a new time slice
            # at the 1/2 steps which is then retimed at the 2/2 steps.
            dd.global_time = time0 # this is the --end-- time, the one we are working on
            substep(actor, Val(:time_advance), δt / 2; progr, retime_core_profiles=(phase == 2))

            if phase == 1
                substep(actor, Val(:evolve_j_ohmic), kk == 1 ? δt / 2 : δt; progr)

                substep(actor, Val(:run_hcd), kk == 1 ? δt / 2 : δt; progr)
            else
                substep(actor, Val(:run_pedestal), kk == 1 ? δt / 2 : δt; progr)

                substep(actor, Val(:run_transport), kk == 1 ? δt / 2 : δt; progr)
                IMAS.time_derivative_source!(dd)
            end

            substep(actor, Val(:run_sawteeth), δt / 2; progr)

            substep(actor, Val(:run_equilibrium), δt / 2; progr)

            substep(actor, Val(:run_pf_active), δt / 2; progr)
        end

    catch e
        @warn "ActorDynamicPlasma failed at $(dd.global_time) [s], trimming the last time slice which might be in a corrupted state"
        IMAS.trim_time!(dd, (-Inf, dd.global_time - 1E-6))
        rethrow(e)
    finally
        actor_logging(dd, old_logging)
    end

    return actor
end

function substep(
    actor::ActorDynamicPlasma,
    ::Val{:time_advance},
    δt::Float64;
    progr=nothing,
    retime_equilibrium::Bool=false,
    retime_core_sources::Bool=false,
    retime_core_profiles::Bool=false
)
    # time_advance
    dd = actor.dd
    if retime_equilibrium
        IMAS.retime!(dd.equilibrium)
    else
        IMAS.new_timeslice!(dd.equilibrium)
    end
    if retime_core_sources
        IMAS.retime!(dd.core_sources)
    else
        IMAS.new_timeslice!(dd.core_sources)
    end
    if retime_core_profiles
        IMAS.retime!(dd.core_profiles)
    else
        IMAS.new_timeslice!(dd.core_profiles)
    end
end

function substep(actor::ActorDynamicPlasma, ::Val{:evolve_j_ohmic}, δt::Float64; progr=nothing, kw...)
    par = actor.par
    # evolve j_ohmic
    if progr !== nothing
        prog, t0, t1, phase = progr
        ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_jt, phase))
    end
    if par.evolve_current
        if actor.actor_jt.par.model == :QED
            if par.ip_controller
                # take vloop from controller
                actor.actor_jt.jt_actor.par.solve_for = :vloop
                actor.actor_jt.jt_actor.par.vloop_from = :controllers__ip
            else
                actor.actor_jt.jt_actor.par.solve_for = :ip
                actor.actor_jt.jt_actor.par.ip_from = :pulse_schedule
            end
            actor.actor_jt.jt_actor.par.Δt = δt
        end
        finalize(step(actor.actor_jt))
    end
end

function substep(actor::ActorDynamicPlasma, ::Val{:run_pedestal}, δt::Float64; progr=nothing, kw...)
    par = actor.par
    # run pedestal actor
    if progr !== nothing
        prog, t0, t1, phase = progr
        ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_ped, phase))
    end
    if par.evolve_pedestal
        finalize(step(actor.actor_ped))
    end
end

function substep(actor::ActorDynamicPlasma, ::Val{:run_transport}, δt::Float64; progr=nothing, kw...)
    par = actor.par
    # run transport actor
    if progr !== nothing
        prog, t0, t1, phase = progr
        ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_tr, phase))
    end
    if par.evolve_transport
        if actor.actor_tr.par.model == :FluxMatcher
            if par.time_derivatives_sources
                actor.actor_tr.tr_actor.par.Δt = δt
            else
                actor.actor_tr.tr_actor.par.Δt = Inf
            end
        end
        finalize(step(actor.actor_tr))
    end
end

function substep(actor::ActorDynamicPlasma, ::Val{:run_equilibrium}, δt::Float64; progr=nothing, kw...)
    par = actor.par
    # run equilibrium actor
    if progr !== nothing
        prog, t0, t1, phase = progr
        ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_eq, phase))
    end
    if par.evolve_equilibrium
        finalize(step(actor.actor_eq))
    end
end

function substep(actor::ActorDynamicPlasma, ::Val{:run_hcd}, δt::Float64; progr=nothing, kw...)
    par = actor.par
    # run HCD actor
    if progr !== nothing
        prog, t0, t1, phase = progr
        ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_hc, phase))
    end
    if par.evolve_hcd
        finalize(step(actor.actor_hc))
    end
end

function substep(actor::ActorDynamicPlasma, ::Val{:run_sawteeth}, δt::Float64; progr=nothing, kw...)
    # run sawteeth actor
    if progr !== nothing
        prog, t0, t1, phase = progr
        ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_saw, phase))
    end
    return finalize(step(actor.actor_saw))
end

function substep(actor::ActorDynamicPlasma, ::Val{:run_pf_active}, δt::Float64; progr=nothing, kw...)
    par = actor.par
    # run pf_active actor
    if progr !== nothing
        prog, t0, t1, phase = progr
        ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_pf, phase))
    end
    if par.evolve_pf_active
        finalize(step(actor.actor_pf))
    end
end

function progress_ActorDynamicPlasma(t0::Float64, t1::Float64, actor::AbstractActor, phase::Int)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]
    return (
        ("    start time", t0),
        ("      end time", t1),
        ("          time", dd.global_time),
        ("         stage", "$(name(actor)) ($phase/2)"),
        ("       Ip [MA]", IMAS.get_from(dd, Val(:ip), :core_profiles) / 1E6),
        ("     Ti0 [keV]", cp1d.t_i_average[1] / 1E3),
        ("     Te0 [keV]", cp1d.electrons.temperature[1] / 1E3),
        ("ne0 [10²⁰ m⁻³]", cp1d.electrons.density_thermal[1] / 1E20),
        ("     max(zeff)", maximum(cp1d.zeff)),
        ("   ω0 [krad/s]", cp1d.rotation_frequency_tor_sonic[1] / 1E3)
    )
end

"""
    plot_plasma_overview(dd::IMAS.dd, time0::Float64)

Plot layout useful to get an overview of a (time dependent) plasma simulation.

Animation in time gets typically done this way:

    #a = @animate for k in 1:length(dd.core_profiles.time)
    @manipulate for k in 1:length(dd.core_profiles.time)
        time0 = dd.core_profiles.time[k]
        FUSE.plot_plasma_overview(dd, time0)
        #savefig("frame_\$(lpad(k-1, 4, '0')).png")
    end
    g = gif(a, "ITER_time_dep.gif", fps=24)
    display(g)

Inclusinon in BEAMER presentation can then be done with:

    \\animategraphics[loop,autoplay,controls,poster=0,width=\\linewidth]{24}{frame_}{0000}{0120}
"""
function plot_plasma_overview(dd::IMAS.dd, time0::Float64=dd.global_time;
    min_power::Float64=0.0,
    aggregate_radiation::Bool=true,
    aggregate_hcd::Bool=true,
    rotation_quantity::Symbol=:sonic,
    dd1::Union{Nothing,IMAS.DD}=nothing,
    size=(1800, 1000),
    kw...)

    GC.safepoint()
    old_gc = GC.enable(false)

    try
        cp1d = dd.core_profiles.profiles_1d[time0]

        # equilibrium cx
        plot_eq_cx = plot()
        if !isempty(dd.build.layer)
            plot!(dd.build; time0, legend=false, equilibrium=false, pf_active=false, label="")
        end
        if dd1 !== nothing
            plot!(dd1.equilibrium.time_slice[time0]; cx=true, color=:black)
        end
        if dd !== dd1
            plot!(dd.equilibrium.time_slice[time0]; cx=true)
        end
        plot!(dd.pf_active; time0, colorbar=nothing)
        plot!(dd.pulse_schedule.position_control, nothing; time0, color=:red)
        out = convex_outline(dd.pf_active.coil)
        if !isempty(out.r)
            plot!(; xlim=[0.0, maximum(out.r)], ylim=extrema(out.z))
        end
        for cw in dd.waves.coherent_wave
            if !isempty(cw.beam_tracing) && time0 >= cw.beam_tracing[1].time
                plot!(cw.beam_tracing[time0]; time0, alpha=0.5, label="")
            end
        end
        plot!(; title="Time = $(@sprintf("%.4f", time0)) [s]", ylabel="", xlabel="")

        # equilibrium top view
        plot_eq_top = plot(; framestyle=:origin, guidefontalign=:bottom)
        plot!(dd.equilibrium.time_slice[time0]; cx=true, top=true)
        for cw in dd.waves.coherent_wave
            if !isempty(cw.beam_tracing) && time0 >= cw.beam_tracing[1].time
                plot!(cw.beam_tracing[time0]; alpha=0.5, top=true)
            end
        end
        if !isempty(dd.wall)
            plot!(dd.wall; top=true)
        end
        plot!(; ylabel="", xlabel="")

        # Ip and Vloop
        plot_ip = plot(; title="Ip [MA] - Vloop [V]")
        if dd1 === nothing
            if !ismissing(dd.pulse_schedule.flux_control, :time)
                plot!(dd.pulse_schedule.flux_control.time,
                    dd.pulse_schedule.flux_control.i_plasma.reference / 1E6;
                    seriestype=:time,
                    color=:black,
                    label="Ip reference [MA]",
                    lw=2.0,
                    ls=:dash,
                    legend_position=:left,
                    background_color_legend=Plots.Colors.RGBA(1.0, 1.0, 1.0, 0.6),
                    legend_foreground_color=:transparent
                )
            end
        else
            plot!(
                dd1.core_profiles.time,
                dd1.core_profiles.global_quantities.ip / 1E6;
                ylim=extrema(dd1.core_profiles.global_quantities.ip / 1E6),
                seriestype=:time,
                color=:black,
                primary=false
            )
        end
        plot!(
            dd.core_profiles.time,
            dd.core_profiles.global_quantities.ip / 1E6;
            ylim=extrema(dd.core_profiles.global_quantities.ip / 1E6),
            seriestype=:time,
            color=:blue,
            label="Ip  [MA]",
            lw=2.0
        )
        # plot!(
        #     dd.equilibrium.time,
        #     [eqt.global_quantities.ip for eqt in dd.equilibrium.time_slice] / 1E6;
        #     ylim=extrema([eqt.global_quantities.ip for eqt in dd.equilibrium.time_slice]/ 1E6),
        #     seriestype=:time,
        #     color=:magenta,
        #     label="Ip  [MA]",
        #     lw=2.0
        # )
        plot!([NaN], [NaN]; color=:red, label="Vloop [V]", lw=2.0)
        vline!([time0]; label="", color=:gray)
        #plot_vloop = plot!(twinx())[2]

        tx = twinx()
        # if dd1 !== nothing
        #     time, data = IMAS.vloop_time(dd1.core_profiles, dd1.equilibrium; method=:area)
        #     plot!(
        #         tx,
        #         time,
        #         data;
        #         seriestype=:time,
        #         color=:black,
        #         label="",
        #         lw=2.0
        #     )
        # end
        if IMAS.controller(dd.controllers, "ip") !== nothing
            time, data = IMAS.vloop_time(dd.controllers)
        else
            time, data = IMAS.vloop_time(dd.core_profiles, dd.equilibrium; method=:area)
        end
        plot!(
            tx,
            time,
            data;
            seriestype=:time,
            color=:red,
            label="",
            lw=2.0
        )

        # core_profiles electrons temperatures
        plot_Te = plot()
        if dd1 !== nothing
            plot!(dd1.core_profiles.profiles_1d[time0].electrons, :temperature; color=:black, only=1, normalization=1E-3, xlabel="", ylabel="", label="Experiment")
            if IMAS.hasdata(dd1.thomson_scattering)
                plot!(dd1.thomson_scattering, :t_e; time0, lw=2.0, normalization=1E-3, xlabel="", ylabel="", label="", primary=false)
            end
        end
        if dd !== dd1
            plot!(cp1d.electrons, :temperature; only=1, lw=2.0, normalization=1E-3, xlabel="", ylabel="", label="Simulated")
        end
        plot!(; title="Electrons temperature [KeV]", xlabel="", ylabel="")

        # core_profiles ions temperatures
        plot_Ti = plot()
        if dd1 !== nothing
            plot!(dd1.core_profiles.profiles_1d[time0], :t_i_average; color=:black, only=1, normalization=1E-3, xlabel="", ylabel="", label="")
            if IMAS.hasdata(dd1.charge_exchange)
                plot!(dd1.charge_exchange, :t_i; time0, lw=2.0, normalization=1E-3, xlabel="", ylabel="", label="", primary=false)
            end
        end
        if dd !== dd1
            plot!(cp1d, :t_i_average; only=1, lw=2.0, normalization=1E-3, xlabel="", ylabel="", label="")
        end
        plot!(; title="Ions temperature [KeV]", xlabel="", ylabel="", label="")

        # core_profiles densities
        plot_n = plot()
        if dd1 !== nothing
            plot!(dd1.core_profiles.profiles_1d[time0]; color=:black, only=2)
            if IMAS.hasdata(dd1.thomson_scattering)
                plot!(dd1.thomson_scattering, :n_e; time0, lw=2.0, xlabel="", ylabel="", label="", primary=false)
            end
            if IMAS.hasdata(dd1.charge_exchange)
                plot!(
                    dd1.charge_exchange,
                    :n_imp;
                    time0,
                    lw=2.0,
                    xlabel="",
                    ylabel="",
                    label="",
                    normalization=dd1.core_profiles.profiles_1d[time0].ion[2].element[1].z_n,
                    primary=false
                )
            end
        end
        if dd !== dd1
            plot!(cp1d; only=2, lw=2.0, legend=:left)
        end
        plot!(; title="Densities [m⁻³]", xlabel="", ylabel="", label="", legend_foreground_color=:transparent)

        # core_profiles rotation
        plot_rotation = plot()
        @assert rotation_quantity in (:toroidal, :sonic)
        if rotation_quantity == :sonic
            if dd1 !== nothing
                if !ismissing(dd1.core_profiles.profiles_1d[time0], :rotation_frequency_tor_sonic)
                    plot!(dd1.core_profiles.profiles_1d[time0], :rotation_frequency_tor_sonic; color=:black, only=1, normalization=1E-3, xlabel="", ylabel="", label="")
                end
                if IMAS.hasdata(dd1.charge_exchange)
                    plot!(dd1.charge_exchange, :ω_tor; time0, lw=2.0, normalization=1E-3, xlabel="", ylabel="", label="", primary=false)
                end
            end
            if dd !== dd1
                if !ismissing(cp1d, :rotation_frequency_tor_sonic)
                    plot!(cp1d, :rotation_frequency_tor_sonic; only=1, lw=2.0, normalization=1E-3, xlabel="", ylabel="", label="")
                end
            end
        else
            if dd1 !== nothing
                if !ismissing(dd1.core_profiles.profiles_1d[time0].ion[1], :rotation_frequency_tor)
                    plot!(dd1.core_profiles.profiles_1d[time0].ion[1], :rotation_frequency_tor; color=:black, only=1, normalization=1E-3, xlabel="", ylabel="", label="")
                end
                if IMAS.hasdata(dd1.charge_exchange)
                    plot!(dd1.charge_exchange, :ω_tor; time0, lw=2.0, normalization=1E-3, xlabel="", ylabel="", label="", primary=false)
                end
            end
            if dd !== dd1
                if !ismissing(cp1d.ion[1], :rotation_frequency_tor)
                    plot!(cp1d.ion[1], :rotation_frequency_tor; only=1, lw=2.0, normalization=1E-3, xlabel="", ylabel="", label="")
                end
            end
        end
        plot!(; title="Toroidal Rotation [Krad/s]", xlabel="", ylabel="", label="")

        # ========

        # core_sources
        plot_j = plot()
        plot!(
            dd.core_sources;
            time0,
            only=5,
            min_power,
            aggregate_radiation,
            aggregate_hcd,
            legend=:topleft,
            legend_foreground_color=:transparent,
            title="Parallel current source [MA/m²]",
            normalization=1E-6,
            ylabel="",
            xlabel=""
        )

        plot_Qe = plot()
        plot!(
            dd.core_sources;
            time0,
            only=1,
            min_power,
            aggregate_radiation,
            aggregate_hcd,
            legend=:topleft,
            legend_foreground_color=:transparent,
            title="Electrons power source [MW/m³]",
            normalization=1E-6,
            ylabel="",
            xlabel=""
        )

        plot_Qi = plot()
        plot!(
            dd.core_sources;
            time0,
            only=2,
            min_power,
            aggregate_radiation,
            aggregate_hcd,
            legend=:topleft,
            legend_foreground_color=:transparent,
            title="Ions power source [MW/m³]",
            normalization=1E-6,
            ylabel="",
            xlabel=""
        )

        plot_Ga = plot()
        plot!(
            dd.core_sources;
            time0,
            only=3,
            min_power,
            aggregate_radiation,
            aggregate_hcd,
            legend=:topleft,
            legend_foreground_color=:transparent,
            title="Particle source [s⁻¹/m³]",
            ylabel="",
            xlabel=""
        )

        plot_Pi = plot()
        plot!(
            dd.core_sources;
            time0,
            only=4,
            min_power,
            aggregate_radiation,
            aggregate_hcd,
            legend=:topleft,
            legend_foreground_color=:transparent,
            title="Rotation sources [Nm]",
            ylabel="",
            xlabel=""
        )

        # ========

        # q
        plot_q = plot(; title="Safety factor q")
        if dd1 !== nothing
            plot!(dd1.equilibrium.time_slice[time0].profiles_1d, :q; color=:black, coordinate=:rho_tor_norm, label="Experiment q")
        end
        sinq = sign(sum(dd.equilibrium.time_slice[time0].profiles_1d.q))
        plot!(dd.equilibrium.time_slice[time0].profiles_1d, :q; lw=2.0, coordinate=:rho_tor_norm, label="Modeled q")
        hline!([sinq * 1]; label="", ls=:dash, color=:black)
        if sinq > 0
            plot!(; ylim=(0, 5))
        else
            plot!(; ylim=(-5, 0))
        end

        # hcd
        plot_hcd = plot(; title="Injected power [MW]")
        plot!(dd.summary.heating_current_drive; ylabel="")
        vline!([time0]; label="", color=:gray)

        # # fusion power
        plot_fusion = plot()
        # plot!(dd.summary.fusion.power, :value; lw=2.0, label="", ylabel="Fusion power [MW]", normalization=5*1E-6, title="Fusion power")
        # vline!([time0]; label="")


        # transport
        plot_trQe = plot()
        plot!(dd.core_transport; time0, only=1, legend=:topleft, legend_foreground_color=:transparent)
        plot_trQi = plot()
        plot!(dd.core_transport; time0, only=2, legend=:topleft, legend_foreground_color=:transparent)
        plot_trGa = plot()
        plot!(dd.core_transport; time0, only=3, legend=:topleft, legend_foreground_color=:transparent)
        plot_trPi = plot()
        plot!(dd.core_transport; time0, only=4, legend=:topleft, legend_foreground_color=:transparent)

        # # inverse scale lengths
        # max_scale = 5
        # subplot = 14
        # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature, :backward); lw=2.0)
        # subplot = 15
        # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.t_i_average, :backward); ylim=(-max_scale, max_scale), lw=2.0)
        # subplot = 16
        # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, :backward); ylim=(-max_scale, max_scale), lw=2.0)

        l = @layout grid(3, 6)
        kw = Dict(kw...)
        mm = Plots.Measures.mm
        p = plot(
            plot_eq_cx, plot_ip, plot_Te, plot_Ti, plot_n, plot_rotation,
            plot_eq_top, plot_j, plot_Qe, plot_Qi, plot_Ga, plot_Pi,
            plot_hcd, plot_q, plot_trQe, plot_trQi, plot_trGa, plot_trPi,
            ; layout=l, left_margin=0 * mm, bottom_margin=0 * mm, size, kw...)

        return p
    finally
        GC.enable(old_gc)
    end
end

function plot_summary(dd::IMAS.DD; kw...)
    return plot_summary(dd.summary; kw...)
end

function plot_summary(summary::IMAS.summary; kw...)
    for leaf in AbstractTrees.Leaves(summary)
        if typeof(leaf.value) <: Vector
            display(plot(leaf.ids, leaf.field; title=IMAS.location(leaf.ids, leaf.field), kw...))
        end
    end
end