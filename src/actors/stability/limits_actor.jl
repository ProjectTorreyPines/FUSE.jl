#= =========== =#
#  ActorStabilityLimits  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorStabilityLimits{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorStabilityLimits<: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorStabilityLimits
    beta_actor::PlasmaAbstractActor
    current_actor::PlasmaAbstractActor
    density_actor::PlasmaAbstractActor
end

"""
ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs all the limit actors. 
"""
function ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorStabilityLimits(kw...)
    actor = ActorStabilityLimits(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorStabilityLimits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, act::ParametersAllActors; kw...)
    logging_actor_init(ActorStabilityLimits)
    par = par(kw...)

    beta_actor = ActorBetaLimit(dd, act)
    current_actor = ActorCurrentLimit(dd, act)
    density_actor = ActorDensityLimit(dd, act)

    return ActorStabilityLimits(dd, par, beta_actor, current_actor, density_actor)
end


"""
    step(actor::ActorStabilityLimits)

Runs through the selected stability actor's step
"""
function _step(actor::ActorStabilityLimits)
    step(actor.beta_actor)
    step(actor.current_actor)
    step(actor.density_actor)
    return actor
end

"""
    finalize(actor::ActorStabilityLimits)

    Finalizes the selected stability actor
"""
function _finalize(actor::ActorStabilityLimits)
    finalize(actor.beta_actor)
    finalize(actor.current_actor)
    finalize(actor.density_actor)
    return actor
end

