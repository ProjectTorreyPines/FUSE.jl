#= ==================== =#
#  ActorStabilityLimits  #
#= ==================== =#
Base.@kwdef mutable struct FUSEparameters__ActorStabilityLimits{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    models::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "Models used for checking plasma stability limits: $(supported_stability_models())"; default=[:default_limits])
    raise_on_breach::Entry{Bool} = Entry{Bool}("-", "Raise an error when one or more stability limits are breached"; default=true)
end

mutable struct ActorStabilityLimits{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorStabilityLimits{P}
    function ActorStabilityLimits(dd::IMAS.dd{D}, par::FUSEparameters__ActorStabilityLimits{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorStabilityLimits)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs all the limit actors
"""
function ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorStabilityLimits(dd, act.ActorStabilityLimits; kw...)
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

    if !isempty(par.models) && par.raise_on_breach
        failed = String[]
        desc = String[]
        for model in dd.stability.model
            if !Bool(@ddtime model.cleared)
                model_name = IMAS.index_2_name__stability__model[model.identifier.index]
                push!(failed, "$(model_name)")
                push!(desc, "$(model_name) ($(@ddtime(model.fraction)) of limit): $(model.identifier.description)")
            end
        end
        if !isempty(failed)
            error("Some stability models have breached their limit threshold:\n$(join(failed, " "))\n* $(join(desc, "\n* "))")
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
