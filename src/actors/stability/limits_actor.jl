#= ==================== =#
#  ActorStabilityLimits  #
#= ==================== =#
Base.@kwdef mutable struct FUSEparameters__ActorStabilityLimits{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    models::Entry{Vector{Symbol}} = Entry(Vector{Symbol}, "-", "Models used for checking plasma stability limits: [$(join(values(IMAS.index_2_name__stability__collection),", ")), $(join(values(IMAS.index_2_name__stability__model),", "))]"; default=[:default_limits])
    raise_on_breach::Entry{Bool} = Entry(Bool, "-", "Raise an error when one or more stability limits are breached"; default=true)
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

Runs all the limit actors
"""
function ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorStabilityLimits
    actor = ActorStabilityLimits(dd, par; kw...)
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
    
    # run all stability models
    run_stability_models(dd, par.models)

    if par.raise_on_breach
        failed = String[]
        time_index = findfirst(dd.stability.time .== @ddtime(dd.stability.time))
        for model in dd.stability.model
            if !Bool(model.cleared[time_index])
                model_name = IMAS.index_2_name__stability__model[model.identifier.index]
                push!(failed, "$(model_name): $(model.identifier.description)")
            end
        end
        if !isempty(failed)
            error("Some stability models have breached their limit threshold:\n* $(join(failed, "\n* "))")
        end
    end

    return actor
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
