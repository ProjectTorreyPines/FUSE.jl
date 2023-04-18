#= ==================== =#
#  ActorStabilityLimits  #
#= ==================== =#

Base.@kwdef mutable struct FUSEparameters__ActorStabilityLimits{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
   model_ids::Entry{Vector{Symbol}} = Entry(Vector{Symbol}, "-", "Models for the limit calculation"; default=[:default])
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
    par = act.ActorStabilityLimits(kw...)
    actor = ActorStabilityLimits(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

"""
    step(actor::ActorStabilityLimits)

Runs through the selected stability actor's step
"""
function _step(actor::ActorStabilityLimits)
    dd = actor.dd
    par = actor.par

    # IS `item` CORRECT? WHAT TO CALL IT?
    for item_id in par.model_ids
        item_index = IMAS.name_2_index(dd.stability.model)[item_id]
        if item_index < 100
            item = resize!(dd.stability.collection, "identifier.index" => item_index)
        else
            item = resize!(dd.stability.model, "identifier.index" => item_index)
        end
        limit_models[item_index](dd, par, item)
    end    

    return actor
end

"""
    finalize(actor::ActorStabilityLimits)

    Finalizes the selected stability actor
"""
function _finalize(actor::ActorStabilityLimits)
    sort!(actor.dd.stability.collection, by=x -> x.identifier.index)
    sort!(actor.dd.stability.model, by=x -> x.identifier.index)
    #actor.dd.stability.model[].cleared
    return actor
end





#######################################################


