#= =================== =#
#  ActorStability       #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorStability{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    stability_actor::Switch{Symbol} = Switch(Symbol, [:BetaLimit, :CurrentLimit, :DensityLimit, :Limits, :None], "-", "Stability Actor to run"; default=:Limits)
    stability_value::Entry{T} = Entry(T, "-", "value of the stability metric"; default=0.0)
    #stability_pass::Entry{Bool} = Entry(Bool, "-", "True if all stability metrics are met"; default=true) #don't really like this, but might be useful
end

mutable struct ActorStability<: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorStability
    stab_actor::PlasmaAbstractActor
    #stab_actor::Union{ActorBetaLimit, ActorCurrentLimit, ActorDensityLimit}
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

    if par.stability_actor == :Limits
        #act.ActorBetaLimit.limit = par.stability_limit
        stab_actor = ActorBetaLimit(dd, act.ActorBetaLimit)
    elseif par.stability_actor == :BetaLimit
        #act.ActorBetaLimit.limit = par.stability_limit 
        stab_actor = ActorBetaLimit(dd, act.ActorBetaLimit)
    else
        error("stability_actor $(par.stability_actor) is not supported yet")
    end

    #value = stab_actor.value
    #add_limit_check!(exceeded_limits, par.stability_actor, value, >, limit)

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