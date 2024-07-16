using ProgressMeter: ProgressMeter
ProgressMeter.ijulia_behavior(:clear)

#= ================== =#
#  ActorDynamicPlasma  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorDynamicPlasma{T<:Real} <: ParametersActorPlasma{T}
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
    #== display and debugging parameters ==#
    verbose::Entry{Bool} = Entry{Bool}("-", "Verbose"; default=false)
end

mutable struct ActorDynamicPlasma{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorDynamicPlasma{P}
    act::ParametersAllActors
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

    if act.ActorCoreTransport.model == :FluxMatcher
        actor_ped = ActorPedestal(
            dd,
            act.ActorPedestal,
            act;
            ip_from=:core_profiles,
            βn_from=:equilibrium, 
            ne_from=:pulse_schedule,
            zeff_ped_from=:pulse_schedule,
            rho_nml=actor_tr.tr_actor.par.rho_transport[end-1],
            rho_ped=actor_tr.tr_actor.par.rho_transport[end])
    else
        actor_ped = ActorNoOperation(dd, act.ActorNoOperation)
    end

    actor_hc = ActorHCD(dd, act.ActorHCD, act)

    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; model=:QED, ip_from=:pulse_schedule, vloop_from=:pulse_schedule)

    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)

    actor_pf = ActorPFactive(dd, act.ActorPFactive)

    return ActorDynamicPlasma(dd, par, act, actor_tr, actor_ped, actor_hc, actor_jt, actor_eq, actor_pf)
end

function _step(actor::ActorDynamicPlasma)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    # time constants
    δt = par.Δt / par.Nt
    t0 = dd.global_time
    t1 = t0 + par.Δt

    # set Δt of the time-dependent actors
    actor.actor_jt.jt_actor.par.Δt = δt
    if actor.actor_tr.par.model == :FluxMatcher
        actor.actor_tr.tr_actor.par.Δt = δt
    end

    # setup things for Ip control
    if par.ip_controller
        η_avg = trapz(cp1d.grid.area, 1.0 ./ cp1d.conductivity_parallel) / cp1d.grid.area[end]
        ctrl_ip = resize!(dd.controllers.linear_controller, "name" => "ip")
        IMAS.pid_controller(ctrl_ip, η_avg * 5.0, η_avg * 0.5, 0.0)
        if IMAS.fxp_request_service(ctrl_ip)
            @warn("Running ip controller service")
        end
        actor.actor_jt.jt_actor.par.solve_for = :vloop
        actor.actor_jt.jt_actor.par.vloop_from = :controllers__ip
    else
        actor.actor_jt.jt_actor.par.solve_for = :ip
        actor.actor_jt.jt_actor.par.ip_from = :pulse_schedule
    end

    prog = ProgressMeter.Progress(par.Nt * 9; dt=0.0, showspeed=true, enabled=par.verbose)
    old_logging = actor_logging(dd, false)

    # this is to fix an issue where having data at the base layer prevents the evaluation of expressions on the lazy copied IDSs
    empty!(dd.core_profiles.profiles_1d[], :j_tor)

    try
        for (kk, tt) in enumerate(range(t0, t1, 2 * par.Nt + 1)[2:end])
            dd.global_time = tt

            # prepare time dependent arrays of structures
            # NOTE: dd.core_profiles is different because it is updated
            #       by actor_jt at the 1/2 steps, but also (mostly)
            #       by actor_tr and actor_ped at the 1/2 steps.
            #       For dd.core_profiles we thus create a new time slice
            #       at the 1/2 steps which is then retimed at the 2/2 steps.
            IMAS.new_timeslice!(dd.equilibrium, tt)
            IMAS.new_timeslice!(dd.core_sources, tt)

            if mod(kk, 2) == 1
                IMAS.new_timeslice!(dd.core_profiles, tt)

                # evolve j_ohmic
                ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_jt, mod(kk, 2) + 1))
                if par.evolve_current
                    if par.ip_controller
                        controller(dd, Val{:ip})
                    end
                    finalize(step(actor.actor_jt))
                end
            else
                IMAS.IMASDD.retime!(dd.core_profiles, tt)

                # run transport actor
                ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_tr, mod(kk, 2) + 1))
                if par.evolve_transport
                    finalize(step(actor.actor_tr))
                end

                # run pedestal actor
                ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_ped, mod(kk, 2) + 1))
                if par.evolve_pedestal
                    finalize(step(actor.actor_ped))
                end
            end

            # run equilibrium actor with the updated beta
            ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_eq, mod(kk, 2) + 1))
            if par.evolve_equilibrium
                finalize(step(actor.actor_eq))
            end

            # run HCD to get updated current drive
            ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_hc, mod(kk, 2) + 1))
            if par.evolve_hcd
                finalize(step(actor.actor_hc))
            end

            # run the pf_active actor to get update coil currents
            ProgressMeter.next!(prog; showvalues=progress_ActorDynamicPlasma(t0, t1, actor.actor_pf, mod(kk, 2) + 1))
            if par.evolve_pf_active
                finalize(step(actor.actor_pf))
            end
        end
    finally
        actor_logging(dd, old_logging)
    end

    return actor
end

function progress_ActorDynamicPlasma(t0::Float64, t1::Float64, actor::AbstractActor, phase::Int)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]
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

Plot layout useful to get an overview of a time dependent plasma simulation.
This gets typically used this way:

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
    p = plot(; layout=l, size=(1600, 1000), margin = 5 * Plots.Measures.mm, kw...)

    cp1d = dd.core_profiles.profiles_1d[time0]

    # Ip and Vloop
    subplot = 1
    plot!(dd.pulse_schedule.flux_control.time,
        dd.pulse_schedule.flux_control.i_plasma.reference / 1E6;
        seriestype=:time,
        color=:gray,
        label="Ip reference [MA]",
        lw=2.0,
        ls=:dash,
        legend_position=:bottomright,
        subplot
    )
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
        plot!([NaN], [NaN]; seriestype=:time, color=:red, label="Vloop [mV]", subplot)
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
    #plot!(dd.equilibrium.time_slice[2]; cx=true, color=:gray, subplot)
    plot!(dd.equilibrium.time_slice[time0]; cx=true, subplot)
    plot!(dd.build; legend=false, subplot)
    plot!(dd.pf_active; time0, subplot, colorbar=:none, axis=false)

    # core_profiles temperatures
    subplot = 3
    plot!(dd.core_profiles.profiles_1d[1]; only=1, color=:gray, label=" initial", subplot, normalization=1E-3)
    plot!(cp1d; only=1, lw=2.0, subplot, normalization=1E-3, ylabel="[keV]")#, ylim=(0.0, 25.0))

    # core_profiles densities
    subplot = 4
    plot!(dd.core_profiles.profiles_1d[1]; only=2, color=:gray, label=" initial", subplot)
    plot!(cp1d; only=2, lw=2.0, subplot, ylabel="[m⁻³]")#, ylim=(0.0, 1.3E20))

    # q
    subplot = 9
    plot!(dd.equilibrium.time_slice[2].profiles_1d, :q; lw=2.0, coordinate=:rho_tor_norm, label="Initial q", subplot)
    plot!(dd.equilibrium.time_slice[time0].profiles_1d, :q; lw=2.0, coordinate=:rho_tor_norm, label="q", subplot)

    # # fusion power
    # subplot = 9
    # plot!(dd.summary.fusion.power, :value; lw=2.0, label="", ylabel="Fusion power [MW]", normalization=5*1E-6, subplot, title="Fusion power")
    # vline!([time0]; label="", subplot)

    # core_sources
    subplot = 5
    plot!(dd.core_sources; time0, only=4, subplot, min_power, aggregate_radiation, weighted=:area, title="Parallel current source", normalization=1E-6, ylabel="[MA]")#, ylim=(0.0, 10.0))

    subplot = 6
    plot!(        
        dd.core_sources;
        time0,
        only=1,
        subplot,
        min_power,
        aggregate_radiation,
        weighted=:volume,
        legend=:bottomleft,
        title="Electron power source",
        normalization=1E-6,
        ylabel="[MW]"
    )#, ylim=(-40.0, 41.0))

    subplot = 7
    plot!(
        dd.core_sources;
        time0,
        only=2,
        subplot,
        min_power,
        aggregate_radiation,
        weighted=:volume,
        legend=:bottomleft,
        title="Ion power source",
        normalization=1E-6,
        ylabel="[MW]"
    )#, ylim=(-40.0, 41.0))

    subplot = 8
    plot!(dd.core_sources; time0, only=3, subplot, min_power, aggregate_radiation, weighted=:volume, title="Electron particle source", ylabel="[s⁻¹]")#, ylim=(0.0, 1.1E20))

    # transport
    #subplot=9
    #plot!(dd.core_transport; time0, only=4, subplot)
    subplot = 10
    plot!(dd.core_transport; time0, only=1, subplot)#, ylim=(0.0, 2.2E5))
    subplot = 11
    plot!(dd.core_transport; time0, only=2, subplot)#, ylim=(0.0, 2.2E5))
    subplot = 12
    plot!(dd.core_transport; time0, only=3, subplot)#, ylim=(0.0, 6.5E17))

    # # inverse scale lenghts
    # max_scale = 5
    # subplot = 14
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature, :backward); subplot, ylim=(-max_scale, max_scale), lw=2.0)
    # subplot = 15
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.t_i_average, :backward); subplot, ylim=(-max_scale, max_scale), lw=2.0)
    # subplot = 16
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, :backward); subplot, ylim=(-max_scale, max_scale), lw=2.0)

    return p
end
