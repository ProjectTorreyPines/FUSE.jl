#= ================= =#
#  ActorCurrentLimit  #
#= ================= =#

Base.@kwdef mutable struct FUSEparameters__ActorCurrentLimit{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:none, :q95], "-", "Model for the current limit calculation"; default=:q95)
end

mutable struct ActorCurrentLimit <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorCurrentLimit
    lim::IMAS.stability__limit
end

"""
    ActorCurrentLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates the emperical current limit driven instabilities using various models
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
    lim = resize!(dd.stability.limit, "name" => "Current limit")
    ActorCurrentLimit(dd, par, lim)
end

"""
    step(actor::ActorCurrentLimit)

Runs ActorCurrentLimit to evaluate the current limit for the given equilibrium
"""
function _step(actor::ActorCurrentLimit)
    dd = actor.dd
    par = actor.par
    lim = actor.lim

    if par.model == :none
        logging(Logging.Debug, :actors, "ActorCurrentLimit: limit check disabled")
    elseif par.model == :standard
        current_q95(dd, par, lim)
    else
        error("ActorCurrentLimit: model = $(par.model) is unknown")
    end

    return actor
end

function _finalize(actor::ActorCurrentLimit)
    lim = actor.lim
    if lim.model.fraction < 1
        lim.cleared = 1
    else
        lim.cleared = 0
    end

    return actor
end

##### Limit Model Functions #####


"""
    current_q95(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__limit)

Standard limit in edge current via the safety factor
Model Formulation: q95 < 2
Citation: 
"""
function current_q95(dd::IMAS.dd, par::FUSEparameters__ActorCurrentLimit, lim::IMAS.stability__limit)
    lim.model.name = "Standard::q95"
    lim.model.formula = "q_95 > 2.0"

    q95 = dd.equilibrium.time_slice[].global_quantities.q_95
    model_value = 1.0 / abs(q95)
    target_value = 0.5

    lim.model.fraction = model_value / target_value
end