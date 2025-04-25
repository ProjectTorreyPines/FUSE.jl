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
    par::OverrideParameters{P,FUSEparameters__ActorDynamicPlasma{P}}
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
        ip_from=:core_profiles,
        βn_from=:core_profiles,
        ne_from=:pulse_schedule,
        zeff_from=:pulse_schedule,
        rho_nml,
        rho_ped)

    actor_hc = ActorHCD(dd, act.ActorHCD, act)

    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; ip_from=:pulse_schedule, vloop_from=:pulse_schedule)

    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)

    actor_pf = ActorPFactive(dd, act.ActorPFactive)

    return ActorDynamicPlasma(dd, par, act, actor_tr, actor_ped, actor_hc, actor_jt, actor_eq, actor_pf)
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
        for (kk, time0) in enumerate(range(t0, t1, 2 * Nt + 1)[2:end]) # NOTE: δt is a full step, some actors are called every 1/2 step
            phase = mod(kk + 1, 2) + 1 # phase can be either 1 or 2, we start with 1
            progr = (prog, t0, t1, phase)

            # Prepare time dependent arrays of structures
            # NOTE: dd.core_profiles is different because it is updated
            # by actor_jt at the 1/2 steps, but also
            # by actor_tr and actor_ped at the 1/2 steps.
            # For dd.core_profiles we thus create a new time slice
            # at the 1/2 steps which is then retimed at the 2/2 steps.
            dd.global_time = time0 # this is the --end-- time, the one we are working on
            substep(actor, Val{:time_advance}, δt / 2; progr, retime_core_profiles=(phase == 2))

            if phase == 1
                substep(actor, Val{:evolve_j_ohmic}, kk == 1 ? δt / 2 : δt; progr)
            else
                substep(actor, Val{:run_pedestal}, kk == 1 ? δt / 2 : δt; progr)

                substep(actor, Val{:run_transport}, kk == 1 ? δt / 2 : δt; progr)
                IMAS.time_derivative_source!(dd)
            end

            substep(actor, Val{:run_equilibrium}, δt / 2; progr)

            substep(actor, Val{:run_hcd}, δt / 2; progr)

            substep(actor, Val{:run_pf_active}, δt / 2; progr)
        end

    finally
        actor_logging(dd, old_logging)
    end

    return actor
end

function substep(
    actor::ActorDynamicPlasma,
    ::Type{Val{:time_advance}},
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

function substep(actor::ActorDynamicPlasma, ::Type{Val{:evolve_j_ohmic}}, δt::Float64; progr=nothing, kw...)
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

function substep(actor::ActorDynamicPlasma, ::Type{Val{:run_pedestal}}, δt::Float64; progr=nothing, kw...)
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

function substep(actor::ActorDynamicPlasma, ::Type{Val{:run_transport}}, δt::Float64; progr=nothing, kw...)
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

function substep(actor::ActorDynamicPlasma, ::Type{Val{:run_equilibrium}}, δt::Float64; progr=nothing, kw...)
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

function substep(actor::ActorDynamicPlasma, ::Type{Val{:run_hcd}}, δt::Float64; progr=nothing, kw...)
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

function substep(actor::ActorDynamicPlasma, ::Type{Val{:run_pf_active}}, δt::Float64; progr=nothing, kw...)
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
        ("       Ip [MA]", IMAS.get_from(dd, Val{:ip}, :core_profiles) / 1E6),
        ("     Ti0 [keV]", cp1d.t_i_average[1] / 1E3),
        ("     Te0 [keV]", cp1d.electrons.temperature[1] / 1E3),
        ("ne0 [10²⁰ m⁻³]", cp1d.electrons.density_thermal[1] / 1E20),
        ("     max(zeff)", maximum(cp1d.zeff))
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
function plot_plasma_overview(dd::IMAS.dd{T}, time0::Float64=dd.global_time;
    min_power::Float64=0.0,
    aggregate_radiation::Bool=true,
    aggregate_hcd::Bool=true,
    dd1::Union{Nothing,IMAS.DD}=nothing,
    kw...) where {T<:Real}

    l = @layout grid(3, 4)
    kw = Dict(kw...)
    if :size ∉ keys(kw)
        kw[:size] = (1600, 1000)
    end
    p = plot(; layout=l, left_margin=1 * Plots.Measures.mm, kw...)

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
    # plot!(
    #     dd.equilibrium.time,
    #     [eqt.global_quantities.ip for eqt in dd.equilibrium.time_slice] / 1E6;
    #     ylim=extrema([eqt.global_quantities.ip for eqt in dd.equilibrium.time_slice]/ 1E6),
    #     seriestype=:time,
    #     color=:magenta,
    #     label="Ip  [MA]",
    #     ylabel="Ip [MA]",
    #     lw=2.0,
    #     subplot
    # )
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
    else
        plot!([NaN], [NaN]; color=:red, label="Vloop [mV]", lw=2.0, subplot)
        time, data = IMAS.vloop_time(dd.core_profiles, dd.equilibrium)#; method=:edge)
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
        plot!(dd.build; time0, subplot, legend=false, equilibrium=false, pf_active=false)
    end
    if dd1 !== nothing
        plot!(dd1.equilibrium.time_slice[time0]; cx=true, color=:black, subplot)
    end
    if dd !== dd1
        plot!(dd.equilibrium.time_slice[time0]; cx=true, subplot)
    end
    plot!(dd.pf_active; time0, subplot, colorbar=nothing)
    plot!(dd.pulse_schedule.position_control, nothing; time0, subplot, color=:red)
    out = convex_outline(dd.pf_active.coil)
    if !isempty(out.r)
        plot!(; xlim=[0.0, maximum(out.r)], ylim=extrema(out.z), subplot)
    end
    plot!(; title="Time = $(@sprintf("%.4f", time0)) [s]", subplot)

    # core_profiles temperatures
    subplot = 3
    #    plot!(dd.core_profiles.profiles_1d[1]; only=1, color=:gray, label=" initial", subplot, normalization=1E-3)
    if dd1 !== nothing
        plot!(dd1.core_profiles.profiles_1d[time0]; color=:black, only=1, subplot, normalization=1E-3, label="")
    end
    if dd !== dd1
        plot!(cp1d; only=1, lw=2.0, subplot, normalization=1E-3, ylabel="[keV]", legend_foreground_color=:transparent)
    end

    # core_profiles densities
    subplot = 4
    #    plot!(dd.core_profiles.profiles_1d[1]; only=2, color=:gray, label=" initial", subplot)
    if dd1 !== nothing
        plot!(dd1.core_profiles.profiles_1d[time0]; color=:black, only=2, subplot, label="")
    end
    if dd !== dd1
        plot!(cp1d; only=2, lw=2.0, subplot, ylabel="[m⁻³]", legend=:left, legend_foreground_color=:transparent)
    end

    # q
    subplot = 9
    #    plot!(dd.equilibrium.time_slice[2].profiles_1d, :q; lw=2.0, coordinate=:rho_tor_norm, label="Initial q", subplot)

    # plot!(
    #     dd.equilibrium.time_slice[time0].profiles_1d,
    #     :q;
    #     lw=2.0,
    #     coordinate=:rho_tor_norm,
    #     label="q",
    #     subplot,
    #     legend_foreground_color=:transparent,
    #     title="Safety factor",
    #     legend=:bottomleft
    # )

    plot!(dd.nbi; subplot, smooth_beam_tau=0.1, legend_foreground_color=:transparent)
    vline!([time0]; label="", subplot)

    # # fusion power
    # subplot = 9
    # plot!(dd.summary.fusion.power, :value; lw=2.0, label="", ylabel="Fusion power [MW]", normalization=5*1E-6, subplot, title="Fusion power")
    # vline!([time0]; label="", subplot)

    # core_sources
    subplot = 5
    if dd1 !== nothing
        plot!(dd1.equilibrium.time_slice[time0].profiles_1d, :q; color=:black, coordinate=:rho_tor_norm, label="Experiment q", subplot)
    end
    plot!(dd.equilibrium.time_slice[time0].profiles_1d, :q; lw=2.0, coordinate=:rho_tor_norm, label="Modeled q", subplot)
    hline!([-1]; subplot, label="", ls=:dash, color=:black)
    plot!(; ylim=(-5, 0), subplot)

    # plot!(
    #     dd.core_sources;
    #     time0,
    #     only=5,
    #     subplot,
    #     min_power,
    #     aggregate_radiation,
    #     aggregate_hcd,
    #     legend=:topleft,
    #     legend_foreground_color=:transparent,
    #     title="Parallel current",
    #     normalization=1E-6,
    #     ylabel="[MA/m²]"
    # )

    subplot = 6
    plot!(
        dd.core_sources;
        time0,
        only=1,
        subplot,
        min_power,
        aggregate_radiation,
        aggregate_hcd,
        legend=:topleft,
        legend_foreground_color=:transparent,
        title="Electron power source",
        normalization=1E-6,
        ylabel="[MW/m³]"
    )

    subplot = 7
    plot!(
        dd.core_sources;
        time0,
        only=2,
        subplot,
        min_power,
        aggregate_radiation,
        aggregate_hcd,
        legend=:topleft,
        legend_foreground_color=:transparent,
        title="Ion power source",
        normalization=1E-6,
        ylabel="[MW/m³]"
    )

    subplot = 8
    plot!(
        dd.core_sources;
        time0,
        only=3,
        subplot,
        min_power,
        aggregate_radiation,
        aggregate_hcd,
        legend=:topleft,
        legend_foreground_color=:transparent,
        title="Particle source",
        ylabel="[s⁻¹/m³]"
    )

    # transport
    #subplot=9
    #plot!(dd.core_transport; time0, only=4, subplot)
    subplot = 10
    plot!(dd.core_transport; time0, only=1, subplot, legend=:topleft, legend_foreground_color=:transparent)
    subplot = 11
    plot!(dd.core_transport; time0, only=2, subplot, legend=:topleft, legend_foreground_color=:transparent)
    subplot = 12
    plot!(dd.core_transport; time0, only=3, subplot, legend=:topleft, legend_foreground_color=:transparent)

    # # inverse scale lengths
    # max_scale = 5
    # subplot = 14
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature, :third_order); subplot, lw=2.0)
    # subplot = 15
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.t_i_average, :third_order); subplot, ylim=(-max_scale, max_scale), lw=2.0)
    # subplot = 16
    # plot!(cp1d.grid.rho_tor_norm, -IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, :third_order); subplot, ylim=(-max_scale, max_scale), lw=2.0)

    return p
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