#= ================= =#
#  ActorPlasmaLimits  #
#= ================= =#
Base.@kwdef mutable struct FUSEparameters__ActorPlasmaLimits{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    models::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "Models used for checking plasma operational limits: $(supported_limit_models)"; default=default_limit_models)
    raise_on_breach::Entry{Bool} = Entry{Bool}("-", "Raise an error when one or more operational limits are breached"; default=true)
end

mutable struct ActorPlasmaLimits{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPlasmaLimits{P}
    act::ParametersAllActors{P}
    function ActorPlasmaLimits(dd::IMAS.dd{D}, par::FUSEparameters__ActorPlasmaLimits{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorPlasmaLimits)
        par = par(kw...)
        return new{D,P}(dd, par, act)
    end
end

"""
    ActorPlasmaLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs all the limit actors
"""
function ActorPlasmaLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPlasmaLimits(dd, act.ActorPlasmaLimits, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorPlasmaLimits)

Runs through the selected plasma limits
"""
function _step(actor::ActorPlasmaLimits)
    dd = actor.dd
    par = actor.par
    act = actor.act
    
    # run all limit models
    run_limits_models(dd, act, par.models)

    if !isempty(par.models) && par.raise_on_breach
        failed = String[]
        desc = String[]
        for model in dd.limits.model
            if !Bool(@ddtime model.cleared)
                model_name = model.identifier.name
                push!(failed, model_name)
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
    _finalize(actor::ActorPlasmaLimits)

Finalizes the selected plasma limits
"""
function _finalize(actor::ActorPlasmaLimits)
    sort!(actor.dd.limits.model, by=x -> x.identifier.name)
    return actor
end
