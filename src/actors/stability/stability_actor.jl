#= =================== =#
#  ActorStability       #
#= =================== =#
const stability_actors = [:Limits, :None]
Base.@kwdef mutable struct FUSEparameters__ActorStability{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    stability_actor::Switch{Symbol} = Switch(Symbol, stability_actors, "-", "Stability Actor to run"; default=:Limits)
end

mutable struct ActorStability<: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorStability
    stab_actor::PlasmaAbstractActor
end

"""
    ActorStability(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple stability actors
"""
function ActorStability(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorStability(kw...)
    actor = ActorStability(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorStability(dd::IMAS.dd, par::FUSEparameters__ActorStability, act::ParametersAllActors; kw...)
    logging_actor_init(ActorStability)
    par = par(kw...)

    if par.stability_actor == :None 
        error("stability_actor $(par.stability_actor) is not supported yet")
    elseif par.stability_actor == :Limits
        stab_actor = ActorStabilityLimits(dd, act)
    else
        error("stability_actor $(par.stability_actor) is not supported yet")
    end

    return ActorStability(dd, par, stab_actor)
end

"""
    step(actor::ActorStability)

Runs through the selected stability actor's step
"""
function _step(actor::ActorStability)
    step(actor.stab_actor)
    return actor
end

"""
    finalize(actor::ActorStability)

    Finalizes the selected stability actor
"""
function _finalize(actor::ActorStability)
    finalize(actor.stab_actor)
    return actor
end

