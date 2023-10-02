using ProgressMeter: ProgressMeter
ProgressMeter.ijulia_behavior(:clear)

#= ================== =#
#  ActorDynamicPlasma  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorDynamicPlasma{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution")
    evolve_transport::Entry{Bool} = Entry{Bool}("-", "Evolve the transport"; default=true)
    evolve_hcd::Entry{Bool} = Entry{Bool}("-", "Evolve the heating and current drive"; default=true)
    evolve_current::Entry{Bool} = Entry{Bool}("-", "Evolve the plasma current"; default=true)
    evolve_equilibrium::Entry{Bool} = Entry{Bool}("-", "Evolve the equilibrium"; default=true)
    evolve_pf_active::Entry{Bool} = Entry{Bool}("-", "Evolve the PF currents"; default=true)
end

mutable struct ActorDynamicPlasma{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorDynamicPlasma{P}
    act::ParametersAllActors
    actor_tr::ActorCoreTransport{D,P}
    actor_hc::ActorHCD{D,P}
    actor_jt::ActorCurrent{D,P}
    actor_eq::ActorEquilibrium{D,P}
    actor_pf::ActorPFcoilsOpt{D,P}
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

    actor_hc = ActorHCD(dd, act.ActorHCD, act)

    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; model=:QED, ip_from=:pulse_schedule, vloop_from=:pulse_schedule)
    actor_jt.jt_actor.par.solve_for = :vloop

    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)

    actor_pf = ActorPFcoilsOpt(dd, act.ActorPFcoilsOpt; optimization_scheme=:currents)

    return ActorDynamicPlasma(dd, par, act, actor_tr, actor_hc, actor_jt, actor_eq, actor_pf)
end

function _step(actor::ActorDynamicPlasma)
    dd = actor.dd
    par = actor.par

    δt = par.Δt / par.Nt
    t0 = dd.global_time
    t1 = t0 + par.Δt

    # set Δt of the current actor
    actor.actor_jt.jt_actor.par.Δt = δt

    # setup things for Ip control
    ctrl_ip = resize!(dd.controllers.linear_controller, "name" => "ip")
    cp1d = dd.core_profiles.profiles_1d[]
    η_avg = integrate(cp1d.grid.area, 1.0 ./ cp1d.conductivity_parallel) / cp1d.grid.area[end]
    fill!(ctrl_ip, IMAS.controllers__linear_controller(η_avg * 1.0, η_avg * 0.5, 0.0))

    prog = ProgressMeter.Progress(par.Nt * 6; dt=0.0, showspeed=true)
    backup_actor_logging = logging()[:actors]
    logging(; actors=Logging.Error)
    try
        for (kk, tt) in enumerate(LinRange(t0, t1, 2 * par.Nt + 1)[2:end])
            # prepare time dependent arrays of structures
            IMAS.new_timeslice!(dd.equilibrium, tt)
            IMAS.new_timeslice!(dd.core_profiles, tt)
            IMAS.new_timeslice!(dd.core_sources, tt)
            dd.global_time = tt

            if mod(kk, 2) == 0
                # run transport actor
                ProgressMeter.next!(prog; showvalues=showvalues(t0, t1, actor.actor_tr, mod(kk, 2) + 1))
                if par.evolve_transport
                    _finalize(_step(actor.actor_tr))
                end
            else
                # evolve j_ohmic
                ProgressMeter.next!(prog; showvalues=showvalues(t0, t1, actor.actor_jt, mod(kk, 2) + 1))
                controller(dd, ctrl_ip, Val{:ip})
                if par.evolve_current
                    _finalize(_step(actor.actor_jt))
                end
            end

            # run equilibrium actor with the updated beta
            ProgressMeter.next!(prog; showvalues=showvalues(t0, t1, actor.actor_eq, mod(kk, 2) + 1))
            if par.evolve_equilibrium
                _finalize(_step(actor.actor_eq))
            end

            # run HCD to get updated current drive
            ProgressMeter.next!(prog; showvalues=showvalues(t0, t1, actor.actor_hc, mod(kk, 2) + 1))
            if par.evolve_hcd
                _finalize(_step(actor.actor_hc))
            end
        end
    catch e
        rethrow(e)
    finally
        logging(; actors=backup_actor_logging)
    end

    # run the pf_active actor to get update coil currents
    # NOTE: for the time being this actor works on all time slices at once
    if par.evolve_pf_active
        finalize(step(actor.actor_pf))
    end

    return actor
end

function showvalues(t0::Float64, t1::Float64, actor::AbstractActor, phase::Int)
    dd = actor.dd
    return (
        ("start time", t0),
        ("  end time", t1),
        ("      time", dd.global_time),
        ("     stage", "$(name(actor)) ($phase/2)"),
        ("        ip", "$(IMAS.get_from(dd, Val{:ip}, :core_profiles)/1E6) [MA]")
    )
end