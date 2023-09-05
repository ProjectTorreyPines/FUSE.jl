ProgressMeter.ijulia_behavior(:clear)

#= ================== =#
#  ActorDynamicPlasma  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorDynamicPlasma{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution")
end

mutable struct ActorDynamicPlasma{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorDynamicPlasma{P}
    act::ParametersAllActors
    actor_tr::ActorCoreTransport{D,P}
    actor_hc::ActorHCD{D,P}
    actor_jt::ActorCurrent{D,P}
    actor_eq::ActorEquilibrium{D,P}
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
    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; model=:QED)
    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act)
    return ActorDynamicPlasma(dd, par, act, actor_tr, actor_hc, actor_jt, actor_eq)
end

function _step(actor::ActorDynamicPlasma)
    dd = actor.dd
    par = actor.par

    δt = par.Δt / par.Nt
    t0 = dd.global_time
    t1 = t0 + par.Δt

    # set Δt of the current actor
    actor.actor_jt.jt_actor.par.Δt = δt

    prog = ProgressMeter.Progress(par.Nt * 6; dt=0.0, showspeed=true)
    backup_actor_logging = logging()[:actors]
    logging(; actors=Logging.Error)
    try
        for tt in LinRange(t0, t1, par.Nt + 1)[2:end]
            begin # first 1/2 step (transport)
                # prepare time dependent arrays of structures
                IMAS.new_timeslice!(dd.equilibrium, tt - δt / 2.0)
                IMAS.new_timeslice!(dd.core_profiles, tt - δt / 2.0)
                IMAS.new_timeslice!(dd.core_sources, tt - δt / 2.0)
                IMAS.new_timeslice!(dd.core_transport, tt - δt / 2.0)
                dd.global_time = tt - δt / 2.0

                #controller() --> (PS.beta - CP.beta) ==> NBI.power

                # run transport actor
                ProgressMeter.next!(prog; showvalues=[("start time", t0), ("  end time", t1), ("      time", dd.global_time), ("     stage", "$(name(actor.actor_tr)) 1/2")])
                finalize(step(actor.actor_tr))

                # run equilibrium actor with the updated beta
                ProgressMeter.next!(prog; showvalues=[("start time", t0), ("  end time", t1), ("      time", dd.global_time), ("     stage", "$(name(actor.actor_eq)) 1/2")])
                finalize(step(actor.actor_eq))

                # run HCD to get updated current drive
                ProgressMeter.next!(prog; showvalues=[("start time", t0), ("  end time", t1), ("      time", dd.global_time), ("     stage", "$(name(actor.actor_hc)) 1/2")])
                finalize(step(actor.actor_hc))
            end

            begin # second 1/2 step (current)
                # prepare time dependent arrays of structures
                IMAS.new_timeslice!(dd.equilibrium, tt)
                IMAS.new_timeslice!(dd.core_profiles, tt)
                IMAS.new_timeslice!(dd.core_sources, tt)
                dd.global_time = tt

                #controller() --> (PS.ip - EQ.ip) ==> OH.Vloop

                # evolve j_ohmic
                ProgressMeter.next!(prog; showvalues=[("start time", t0), ("  end time", t1), ("      time", dd.global_time), ("     stage", "$(name(actor.actor_jt)) 2/2")])
                finalize(step(actor.actor_jt))

                # run equilibrium actor with the updated beta
                ProgressMeter.next!(prog; showvalues=[("start time", t0), ("  end time", t1), ("      time", dd.global_time), ("     stage", "$(name(actor.actor_eq)) 2/2")])
                finalize(step(actor.actor_eq))

                # run HCD to get updated current drive
                ProgressMeter.next!(prog; showvalues=[("start time", t0), ("  end time", t1), ("      time", dd.global_time), ("     stage", "$(name(actor.actor_hc)) 2/2")])
                finalize(step(actor.actor_hc))
            end
            println(stderr, dd.global_time)
        end
    catch e
        logging(; actors=backup_actor_logging)
        rethrow(e)
    end

    return actor
end
