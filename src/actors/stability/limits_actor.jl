#= ==================== =#
#  ActorStabilityLimits  #
#= ==================== =#

Base.@kwdef mutable struct FUSEparameters__ActorStabilityLimits{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
   model_ids::Entry{Vector{Symbol}} = Entry(Vector{Symbol}, "-", "Models for the limit calculation"; default=[:force_fail])
end

mutable struct ActorStabilityLimits <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorStabilityLimits
end

"""
ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs all the limit actors. 
"""
function ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorStabilityLimits(kw...)
    actor = ActorStabilityLimits(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorStabilityLimits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits; kw...)
    logging_actor_init(ActorStabilityLimits)
    par = par(kw...)
    @ddtime(dd.stability.time = dd.global_time)
    return ActorStabilityLimits(dd, par)
end

"""
    step(actor::ActorStabilityLimits)

Runs through the selected stability actor's step
"""
function _step(actor::ActorStabilityLimits)
    dd = actor.dd
    par = actor.par

    for model_id in par.model_ids
        model_index = IMAS.name_2_index(dd.stability.model)[model_id]
        model = resize!(dd.stability.model, "identifier.index" => model_index)
        limit_models[model_index](dd, par, model)
    end    

    return actor
end

"""
    finalize(actor::ActorStabilityLimits)

    Finalizes the selected stability actor
"""
function _finalize(actor::ActorStabilityLimits)
    # put a sort here by index
    return actor
end





#######################################################


