#= ================== =#
#  ActorDynamicPlasma  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorDynamicPlasma{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution")
    evolve_transport::Entry{Bool} = Entry{Bool}("-", "Evolve the transport"; default=true)
    evolve_pedestal::Entry{Bool} = Entry{Bool}("-", "Evolve the pedestal"; default=true)
    evolve_hcd::Entry{Bool} = Entry{Bool}("-", "Evolve the heating and current drive"; default=true)
    evolve_current::Entry{Bool} = Entry{Bool}("-", "Evolve the plasma current"; default=true)
    evolve_equilibrium::Entry{Bool} = Entry{Bool}("-", "Evolve the equilibrium"; default=true)
    evolve_pf_active::Entry{Bool} = Entry{Bool}("-", "Evolve the PF currents"; default=true)
    ip_controller::Entry{Bool} = Entry{Bool}("-", "Use controller to change v_loop to match desired Ip"; default=false)
    time_derivatives_sources::Entry{Bool} = Entry{Bool}("-", "Include time-derivative sources"; default=false)
    #== display and debugging parameters ==#
    verbose::Entry{Bool} = Entry{Bool}("-", "Verbose"; default=false)
end

mutable struct ActorDynamicPlasma{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorDynamicPlasma{P}
    act::ParametersAllActors{P}
    actor_tr::ActorCoreTransport{D,P}
    actor_ped::Union{ActorPedestal{D,P},ActorNoOperation{D,P}}
    actor_hc::ActorHCD{D,P}
    actor_jt::ActorCurrent{D,P}
    actor_eq::ActorEquilibrium{D,P}
    actor_pf::ActorPFactive{D,P}
end

"""
    ActorDynamicPlasma(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound evolves plasma in time
"""
function ActorDynamicPlasma(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorDynamicPlasma(dd, act.ActorDynamicPlasma, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorDynamicPlasma(dd::IMAS.dd, par::FUSEparameters__ActorDynamicPlasma, act::ParametersAllActors; kw...)
    logging_actor_init(ActorDynamicPlasma)
    par = par(kw...)

    actor_tr = ActorCoreTransport(dd, act.ActorCoreTransport, act)

    if act.ActorCoreTransport.model in (:FluxMatcher, :EPEDProfiles)
        # allows users to hardwire `rho_nml` and `rho_ped`
        if act.ActorCoreTransport.model == :FluxMatcher && ismissing(act.ActorPedestal, :rho_nml)
            rho_nml = actor_tr.tr_actor.par.rho_transport[end-1]
        else
            rho_nml = act.ActorPedestal.rho_nml
        end
        if act.ActorCoreTransport.model == :FluxMatcher && ismissing(act.ActorPedestal, :rho_ped)
            rho_ped = actor_tr.tr_actor.par.rho_transport[end]
        else
            rho_ped = act.ActorPedestal.rho_ped
        end
        actor_ped = ActorPedestal(
            dd,
            act.ActorPedestal,
            act;
            ip_from=:core_profiles,
            βn_from=:core_profiles,
            ne_from=:pulse_schedule,
            zeff_ped_from=:pulse_schedule,
            rho_nml,
            rho_ped)
    else
        actor_ped = ActorNoOperation(dd, act.ActorNoOperation)
    end

    actor_hc = ActorHCD(dd, act.ActorHCD, act)

    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; model=:QED, ip_from=:pulse_schedule, vloop_from=:pulse_schedule)

    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)

    actor_pf = ActorPFactive(dd, act.ActorPFactive)

    return ActorDynamicPlasma(dd, par, act, actor_tr, actor_ped, actor_hc, actor_jt, actor_eq, actor_pf)
end

function _step(actor::ActorDynamicPlasma; n_steps::Int=0)
    dd = actor.dd
    par = actor.par

    Δt = par.Δt
    Nt = par.Nt
    δt = Δt / Nt

    if n_steps != 0
        Nt = n_steps
        Δt = δt * δt
    end

    # time constants
    t0 = dd.global_time
    t1 = t0 + Δt

    # set Δt of the time-dependent actors
    if actor.actor_jt.par.model == :QED
        actor.actor_jt.jt_actor.par.Δt = δt
    end
    if actor.actor_tr.par.model == :FluxMatcher
        if par.time_derivatives_sources
            actor.actor_tr.tr_actor.par.Δt = δt
        else
            actor.actor_tr.tr_actor.par.Δt = Inf
        end
    end

    # setup things for Ip control
    if par.ip_controller
        ctrl_ip = ip_control(dd, δt)
        # take vloop from controller
        actor.actor_jt.jt_actor.par.solve_for = :vloop
        actor.actor_jt.jt_actor.par.vloop_from = :controllers__ip
    else
        actor.actor_jt.jt_actor.par.solve_for = :ip
        actor.actor_jt.jt_actor.par.ip_from = :pulse_schedule
    end

    step_calls_per_2loop = 9
    ProgressMeter.ijulia_behavior(:clear)
    prog = ProgressMeter.Progress(Nt * step_calls_per_2loop; dt=0.0, showspeed=true, enabled=par.verbose)
    old_logging = actor_logging(dd, false)

    # remove time dependent data after global_time
    IMAS.trim_time!(dd, (-Inf, t0); trim_pulse_schedule=false)

    try
        for (kk, tt) in enumerate(range(t0, t1, 2 * Nt + 1)[2:end])
            phase = mod(kk + 1, 2) + 1 # phase can be either 1 or 2
            # logging(Logging.Info, :actors, " "^workflow_depth(actor.dd) * "--------------- $(Int(ceil(kk/2))) $phase/2 @ $tt")

            dd.global_time = tt

            # prepare time dependent arrays of structures
            # NOTE: dd.core_profiles is different because it is updated
            #       by actor_jt at the 1/2 steps, but also (mostly)
            #       by actor_tr and actor_ped at the 1/2 steps.
            #       For dd.core_profiles we thus create a new time slice
            #       at the 1/2 steps which is then retimed at the 2/2 steps.
            IMAS.new_timeslice!(dd.equilibrium, tt)
            IMAS.new_timeslice!(dd.core_sources, tt)
            if phase == 1
                IMAS.new_timeslice!(dd.core_profiles, tt)
            else
                IMAS.retime!(dd.core_profiles, tt)
            end

            if phase == 1
                # evolve j_ohmic
                ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_jt, phase))
                if par.evolve_current
                    if par.ip_controller
                        ip_control(ctrl_ip, dd)
                    end
                    finalize(step(actor.actor_jt))
                end
            else
                # run pedestal actor
                ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_ped, phase))
                if par.evolve_pedestal
                    finalize(step(actor.actor_ped))
                end

                # run transport actor
                ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_tr, phase))
                if par.evolve_transport
                    finalize(step(actor.actor_tr))
                end
            end

            # run equilibrium actor with the updated beta
            ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_eq, phase))
            if par.evolve_equilibrium
                finalize(step(actor.actor_eq))
            end

            # run HCD to get updated current drive
            ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_hc, phase))
            if par.evolve_hcd
                finalize(step(actor.actor_hc))
            end

            # run the pf_active actor to get update coil currents
            ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_pf, phase))
            if par.evolve_pf_active
                finalize(step(actor.actor_pf))
            end

            # allow running `FUSE.step(actor; n_steps=1)`
            if n_steps > 0 && Int(floor(kk / 2)) == n_steps
                break
            end
        end

    finally
        actor_logging(dd, old_logging)
    end

    # logging(Logging.Info, :actors, " "^workflow_depth(actor.dd) * "---------------")

    return actor
end

function progress_ActorDynamicPlasma(t0::Float64, t1::Float64, actor::AbstractActor, phase::Int)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[end]
    return (
        ("    start time", t0),
        ("      end time", t1),
        ("          time", dd.global_time),
        ("         stage", "$(name(actor)) ($phase/2)"),
        ("       Ip [MA]", IMAS.get_from(dd, Val{:ip}, :core_profiles) / 1E6),
        ("     Ti0 [keV]", cp1d.t_i_average[1] / 1E3),
        ("     Te0 [keV]", cp1d.electrons.temperature[1] / 1E3),
        ("ne0 [10²⁰ m⁻³]", cp1d.electrons.density_thermal[1] / 1E20)
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
function plot_plasma_overview(dd::IMAS.dd, time0::Float64=dd.global_time; min_power::Float64=0.0, aggregate_radiation::Bool=true, kw...)
    l = @layout grid(3, 4)
    p = plot(; layout=l, size=(1600, 1000), left_margin=1 * Plots.Measures.mm, kw...)

    cp1d = dd.core_profiles.profiles_1d[time0]

    # Ip and Vloop
    subplot = 1
    if !ismissing(dd.pulse_schedule.flux_control, :time)
        plot!(dd.pulse_schedule.flux_control.time,
            dd.pulse_schedule.flux_control.i_plasma.reference / 1E6;
            seriestype=:time,
            color=:gray,
            label="Ip reference [MA]",
            lw=2.0,
            ls=:dash,
            legend_position=:left,
            background_color_legend=Plots.Colors.RGBA(1.0, 1.0, 1.0, 0.6),
            legend_foreground_color=:transparent,
            subplot
        )
    end
    plot!(
        dd.core_profiles.time,
        dd.core_profiles.global_quantities.ip / 1E6;
        ylim=extrema(dd.core_profiles.global_quantities.ip / 1E6),
        seriestype=:time,
        color=:blue,
        label="Ip  [MA]",
        ylabel="Ip [MA]",
        lw=2.0,
        subplot
    )
    if IMAS.controller(dd.controllers, "ip") !== nothing
        plot!([NaN], [NaN]; color=:red, label="Vloop [mV]", lw=2.0, subplot)
        time, data = IMAS.vloop_time(dd.controllers)
        plot!(
            twinx(),
            time,
            data .* 1E3;
            seriestype=:time,
            color=:red,
            ylabel="Vloop [mV]",
            label="",
            lw=2.0,
            subplot
        )
    end
    vline!([time0]; label="", subplot)

    # equilibrium, build, and pf_active
    subplot = 2
    if !isempty(dd.build.layer)
        plot!(dd.build; time0, subplot, axis=false, legend=false)
    else
        plot!(dd.equilibrium.time_slice[time0]; cx=true, subplot)
    end
    plot!(dd.pulse_schedule.position_control; time0, subplot, color=:red)
    out = convex_outline(dd.pf_active.coil)
    if !isempty(out.r)
        plot!(; xlim=[0.0, maximum(out.r)], ylim=extrema(out.z), subplot)
    end

    # core_profiles temperatures
    subplot = 3
    #    plot!(dd.core_profiles.profiles_1d[1]; only=1, color=:gray, label=" initial", subplot, normalization=1E-3)
    plot!(cp1d; only=1, lw=2.0, subplot, normalization=1E-3, ylabel="[keV]", legend_foreground_color=:transparent)#, ylim=(0.0, 23.0))

    # core_profiles densities
    subplot = 4
    #    plot!(dd.core_profiles.profiles_1d[1]; only=2, color=:gray, label=" initial", subplot)
    plot!(cp1d; only=2, lw=2.0, subplot, ylabel="[m⁻³]", legend=:left, legend_foreground_color=:transparent)#, ylim=(0.0, 1.3E20))

    # q
    subplot = 9
#    plot!(dd.equilibrium.time_slice[2].profiles_1d, :q; lw=2.0, coordinate=:rho_tor_norm, label="Initial q", subplot)
    plot!(
        dd.equilibrium.time_slice[time0].profiles_1d,
        :q;
        lw=2.0,
        coordinate=:rho_tor_norm,
        label="q",
        subplot,
        legend_foreground_color=:transparent,
        title="Safety factor",
        legend=:bottomleft
    )

    # # fusion power
    # subplot = 9
    # plot!(dd.summary.fusion.power, :value; lw=2.0, label="", ylabel="Fusion power [MW]", normalization=5*1E-6, subplot, title="Fusion power")
    # vline!([time0]; label="", subplot)

    # core_sources
    subplot = 5
    plot!(
        dd.core_sources;
        time0,
        only=5,
        subplot,
        min_power,
        aggregate_radiation,
        weighted=:area,
        legend=:topleft,
        legend_foreground_color=:transparent,
        title="Parallel current",
        normalization=1E-6,
        ylabel="[MA]",
        #ylim=(0.0, 10.0)
    )

    subplot = 6
    plot!(
        dd.core_sources;
        time0,
        only=1,
        subplot,
        min_power,
        aggregate_radiation,
        weighted=:volume,
        legend=:topleft,
        legend_foreground_color=:transparent,
        title="Electron power source",
        normalization=1E-6,
        ylabel="[MW]",
        #ylim=(-20.0, 41.0)
    )

    subplot = 7
    plot!(
        dd.core_sources;
        time0,
        only=2,
        subplot,
        min_power,
        aggregate_radiation,
        weighted=:volume,
        legend=:topleft,
        legend_foreground_color=:transparent,
        title="Ion power source",
        normalization=1E-6,
        ylabel="[MW]",
        #ylim=(-20.0, 41.0)
    )

    subplot = 8
    plot!(
        dd.core_sources;
        time0,
        only=3,
        subplot,
        min_power,
        aggregate_radiation,
        weighted=:volume,
        legend=:topleft,
        legend_foreground_color=:transparent,
        title="Particle source",
        ylabel="[s⁻¹]",
        #ylim=(-0.3E20, 1.1E20)
    )

    # transport
    #subplot=9
    #plot!(dd.core_transport; time0, only=4, subplot)
    subplot = 10
    plot!(dd.core_transport; time0, only=1, subplot, legend=:bottomleft, legend_foreground_color=:transparent)#, ylim=(0.0, 0.3))
    subplot = 11
    plot!(dd.core_transport; time0, only=2, subplot, legend=:bottomleft, legend_foreground_color=:transparent)#, ylim=(0.0, 0.3))
    subplot = 12
    plot!(dd.core_transport; time0, only=3, subplot, legend=:bottomleft, legend_foreground_color=:transparent)#, ylim=(0.0, 6.5E17))

    # # inverse scale lenghts
    # max_scale = 5
    # subplot = 14
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature, :third_order); subplot, ylim=(-max_scale, max_scale), lw=2.0)
    # subplot = 15
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.t_i_average, :third_order); subplot, ylim=(-max_scale, max_scale), lw=2.0)
    # subplot = 16
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, :third_order); subplot, ylim=(-max_scale, max_scale), lw=2.0)

    return p
end

function plot_summary(dd::IMAS.DD; kw...)
    plot_summary(dd.summary; kw...)
end

function plot_summary(summary::IMAS.summary; kw...)
    for leaf in AbstractTrees.Leaves(summary)
        if typeof(leaf.value) <: Vector
            display(plot(leaf.ids, leaf.field; title=IMAS.location(leaf.ids,leaf.field), kw...))
        end
    end
end