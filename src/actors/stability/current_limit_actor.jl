#= ================= =#
#  ActorCurrentLimit  #
#= ================= =#

Base.@kwdef mutable struct FUSEparameters__ActorCurrentLimit{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:Standard, :None], "-", "Model for the limit calculation"; default=:None)
    submodel::Switch{Symbol} = Switch(Symbol, [:A, :B, :C, :None], "-", "Submodel for the limit calculation"; default=:None)
end

mutable struct ActorCurrentLimit <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorCurrentLimit
    limit::IMAS.stability__limit
end

"""
    ActorCurrentLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)


"""
function ActorCurrentLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCurrentLimit(kw...)
    actor = ActorCurrentLimit(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCurrentLimit(dd::IMAS.dd, par::FUSEparameters__ActorCurrentLimit; kw...)
    logging_actor_init(ActorCurrentLimit)
    par = par(kw...)
    limit = resize!(dd.stability.limit, "name" => "Current limit")
    ActorCurrentLimit(dd, par, limit)
end

"""
    step(actor::ActorCurrentLimit)

Runs ActorCurrentLimit to evaluate the Current limit for the given equilibrium
"""
function _step(actor::ActorCurrentLimit)
    dd = actor.dd
    par = actor.par
    limit = actor.limit

    eqt = dd.equilibrium.time_slice[]

    if par.model == :None 
        logging(Logging.Error, :actors, "ActorCurrentLimit: limit check disabled")
    elseif par.model == :Standard
        model_value = 1/abs(eqt.global_quantities.q_95)
        target_value = current_standard_1(dd, par, limit)
    else
        error("ActorCurrentLimit: model = $(par.model) is unknown")
    end

    limit.model.fraction = model_value / target_value
    return actor
end

function _finalize(actor::ActorCurrentLimit)
    limit = actor.limit
    if limit.model.fraction  < 1
        limit.cleared = 1
    else
        limit.cleared = 0 
    end

    return actor
end

##### Limit Model Functions #####

function current_standard_1(dd, par, limit)
    limit.model.name = "Standard::q95"
    limit.model.formula = "q_95 > 2.0"
    target_value = limit.model.limit = 0.5
    return target_value
end