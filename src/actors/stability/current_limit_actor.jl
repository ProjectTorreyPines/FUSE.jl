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

    if par.model == :None 
        logging(Logging.Error, :actors, "ActorCurrentLimit: limit check disabled")
    elseif par.model == :Standard
        current_standard_a(dd, par, lim)
    else
        error("ActorCurrentLimit: model = $(par.model) is unknown")
    end

    return actor
end

function _finalize(actor::ActorCurrentLimit)
    lim = actor.lim
    if lim.model.fraction  < 1
        lim.cleared = 1
    else
        lim.cleared = 0 
    end

    return actor
end

##### Limit Model Functions #####

"""
    current_standard(dd::IMAS.dd)

Standard limit in edge current via the safety factor
Model Formulation: q95 < C
Citation: 
"""
function current_standard(dd::IMAS.dd)
    
    q95 = dd.equilibrium.time_slice[].global_quantities.q_95 
    
    model_value = 1/abs(q95)

    return model_value
end

"""
    current_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__limit)

Standard limit in edge current via the safety factor
Model Formulation: q95 < 2
Citation: 
"""
function current_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorCurrentLimit, lim::IMAS.stability__limit)
    lim.model.name = "Standard::q95"
    lim.model.formula = "q_95 > 2.0"

    model_value = current_standard(dd)
    target_value = 0.5

    lim.model.fraction = model_value / target_value
end