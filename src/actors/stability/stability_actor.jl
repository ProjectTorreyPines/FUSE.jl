#= =================== =#
#  ActorStability       #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorStability{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    stability_actor::Switch{Symbol} = Switch(Symbol, [:StabilityLimits, :none], "-", "Stability Actor to run"; default=:StabilityLimits)
end

mutable struct ActorStability <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorStability
    stab_actor::Union{Nothing,ActorStabilityLimits}
end

function ActorStability(dd::IMAS.dd, par::FUSEparameters__ActorStability, act::ParametersAllActors; kw...)
    logging_actor_init(ActorStability)
    par = par(kw...)

    if par.stability_actor == :none
        stab_actor = nothing
    elseif par.stability_actor == :StabilityLimits
        stab_actor = ActorStabilityLimits(dd, act)
    else
        error("stability_actor `$(par.stability_actor)` is not supported yet")
    end

    return ActorStability(dd, par, stab_actor)
end

"""
    ActorStability(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple stability actors
"""
function ActorStability(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorStability
    actor = ActorStability(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorStability)

Runs through the selected stability actor's step
"""
function _step(actor::ActorStability)
    step(actor.stab_actor)
    return actor
end

"""
    _finalize(actor::ActorStability)

Finalizes the selected stability actor
"""
function _finalize(actor::ActorStability)
    finalize(actor.stab_actor)
    actor.dd.stability.all_cleared
    return actor
end
