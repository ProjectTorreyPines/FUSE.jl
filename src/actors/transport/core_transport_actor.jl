#= ================== =#
#  ActorCoreTransport  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorCoreTransport{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch{Symbol}([:Tauenn, :FluxMatcher, :FixedProfiles], "-", "Transport actor to run"; default=:FixedProfiles)
end

mutable struct ActorCoreTransport{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCoreTransport{P}
    tr_actor::Union{ActorFluxMatcher{D,P},ActorTauenn{D,P},ActorFixedProfiles{D,P}}
end

"""
    ActorCoreTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple core transport actors
"""
function ActorCoreTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCoreTransport(dd, act.ActorCoreTransport, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCoreTransport(dd::IMAS.dd, par::FUSEparameters__ActorCoreTransport, act::ParametersAllActors; kw...)
    logging_actor_init(ActorCoreTransport)
    par = par(kw...)
    if par.model == :FluxMatcher
        tr_actor = ActorFluxMatcher(dd, act.ActorFluxMatcher, act)
    elseif par.model == :Tauenn
        tr_actor = ActorTauenn(dd, act.ActorTauenn)
    elseif par.model == :FixedProfiles
        tr_actor = ActorFixedProfiles(dd, act.ActorFixedProfiles, act)
    end
    return ActorCoreTransport(dd, par, tr_actor)
end

"""
    step(actor::ActorCoreTransport)

Runs through the selected core transport actor step
"""
function _step(actor::ActorCoreTransport)
    step(actor.tr_actor)
    return actor
end

"""
    _finalize(actor::ActorCoreTransport)

Finalizes the selected core transport actor finalize
"""
function _finalize(actor::ActorCoreTransport)
    finalize(actor.tr_actor)
    return actor
end