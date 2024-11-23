#= ================== =#
#  ActorCoreTransport  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorCoreTransport{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:FluxMatcher, :BetaMatch, :EPEDProfiles, :none], "-", "Transport actor to run"; default=:FluxMatcher)
    do_plot::Entry{Bool} = act_common_parameters(do_plot=false)
end

mutable struct ActorCoreTransport{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCoreTransport{P}
    tr_actor::Union{ActorFluxMatcher{D,P},ActorEPEDprofiles{D,P},ActorNoOperation{D,P}}
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
        tr_actor = ActorFluxMatcher(dd, act.ActorFluxMatcher, act; par.do_plot)
    elseif par.model == :EPEDProfiles
        tr_actor = ActorEPEDprofiles(dd, act.ActorEPEDprofiles, act)
    elseif par.model == :BetaMatch
        tr_actor = ActorBetaMatch(dd, act.ActorBetaMatch, act)
    elseif par.model == :none
        tr_actor = ActorNoOperation(dd, act.ActorNoOperation)
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