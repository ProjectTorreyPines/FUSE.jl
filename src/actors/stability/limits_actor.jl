#= ==================== =#
#  ActorStabilityLimits  #
#= ==================== =#

Base.@kwdef mutable struct FUSEparameters__ActorStabilityLimits{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    models::Entry{Vector{Symbol}} = Entry(Vector{Symbol}, "-", "Models for the limit calculation"; default=[:default])
end

mutable struct ActorStabilityLimits <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorStabilityLimits
    function ActorStabilityLimits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits; kw...)
        logging_actor_init(ActorStabilityLimits)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs all the limit actors. 
"""
function ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorStabilityLimits
    actor = ActorStabilityLimits(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorStabilityLimits)

Runs through the selected stability actor's step
"""
function _step(actor::ActorStabilityLimits)
    dd = actor.dd
    par = actor.par
    run_stability_models(dd, par.models)
end

"""
    _finalize(actor::ActorStabilityLimits)

Finalizes the selected stability actor
"""
function _finalize(actor::ActorStabilityLimits)
    sort!(actor.dd.stability.collection, by=x -> x.identifier.index)
    sort!(actor.dd.stability.model, by=x -> x.identifier.index)
    return actor
end
