#= ================== =#
#  ActorDynamicPlasma  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorDynamicPlasma{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution", default=100)
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
    actor_jt = ActorCurrent(dd, act.ActorCurrent, act)
    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act)
    return ActorDynamicPlasma(dd, par, act, actor_tr, actor_hc, actor_jt, actor_eq)
end

function _step(actor::ActorDynamicPlasma)
    dd = actor.dd
    par = actor.par

    δt = par.Δt / par.Nt
    t0 = dd.global_time
    t1 = t0 + par.Δt

    while dd.global_time < t1
        begin # first 1/2 step

            # prepare time dependent arrays of structures
            tt = dd.global_time + δt / 2.0
            IMAS.new_timeslice!(dd.equilibrirum, tt)
            IMAS.new_timeslice!(dd.core_profiles, tt)
            IMAS.new_timeslice!(dd.core_source, tt)
            IMAS.new_timeslice!(dd.core_transport, tt)
            dd.global_time = tt

            # run transport actor
            finalize(step(actor.actor_tr))

            # run equilibrium actor with the updated beta
            finalize(step(actor.actor_eq))

            # run HCD to get updated current drive
            finalize(step(actor.actor_hc))
        end

        begin # second 1/2 step

            # prepare time dependent arrays of structures
            tt = dd.global_time + δt / 2.0
            IMAS.new_timeslice!(dd.equilibrirum, tt)
            IMAS.new_timeslice!(dd.core_profiles, tt)
            IMAS.new_timeslice!(dd.core_source, tt)
            IMAS.new_timeslice!(dd.core_transport, tt)
            dd.global_time = tt

            # evolve j_ohmic
            finalize(step(actor.actor_jt))

            # run equilibrium actor with the updated beta
            finalize(step(actor.actor_eq))

            # run HCD to get updated current drive
            finalize(step(actor.actor_hc))
        end
    end

    return actor
end
