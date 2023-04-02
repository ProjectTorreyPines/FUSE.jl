#= =================== =#
#  ActorStability       #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorStability{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    stability_actor::Switch{Symbol} = Switch(Symbol, [:BetaLimit, :CurrentLimit, :DensityLimit, :Limits, :None], "-", "Stability Actor to run"; default=:Limits)
end

mutable struct ActorStability<: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorStability
    stab_actor::PlasmaAbstractActor
    #stab_actor::Union{Nothing, ActorBetaLimit, ActorCurrentLimit, ActorDensityLimit}
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
    elseif par.stability_actor == :BetaLimit
        stab_actor = ActorBetaLimit(dd, act.ActorBetaLimit)
    elseif par.stability_actor == :CurrentLimit
        stab_actor = ActorCurrentLimit(dd, act.ActorCurrentLimit)
    elseif par.stability_actor == :DensityLimit
        stab_actor = ActorDensityLimit(dd, act.ActorDensityLimit)
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
    lim = actor.dd.stability.limit
    for limit in lim
        if Bool(limit.cleared)
            println("$(limit.name) all clear")
        else
            println("$(limit.name) failed: $(trunc(Int64,limit.model.fraction*100))% of limit")
        end
    end
    return actor
end

