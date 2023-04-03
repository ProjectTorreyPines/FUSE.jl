#= ================= =#
#  ActorDensityLimit  #
#= ================= =#

Base.@kwdef mutable struct FUSEparameters__ActorDensityLimit{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:none, :greenwald], "-", "Model for the limit calculation"; default=:greenwald)
end

mutable struct ActorDensityLimit <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorDensityLimit
    lim::IMAS.stability__limit
end

"""
    ActorDensityLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates the emperical density limit using various models
"""
function ActorDensityLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorDensityLimit(kw...)
    actor = ActorDensityLimit(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorDensityLimit(dd::IMAS.dd, par::FUSEparameters__ActorDensityLimit; kw...)
    logging_actor_init(ActorDensityLimit)
    par = par(kw...)
    lim = resize!(dd.stability.limit, "name" => "Density limit")
    ActorDensityLimit(dd, par, lim)
end


"""
    step(actor::ActorDensityLimit)

Runs ActorDensityLimit to evaluate the Density limit for the given equilibrium
"""
function _step(actor::ActorDensityLimit)
    dd = actor.dd
    par = actor.par
    lim = actor.lim

    if par.model == :none
        logging(Logging.Debug, :actors, "ActorDensityLimit: limit check disabled")
    elseif par.model == :greenwald
        density_greenwald(dd, par, lim)
    else
        error("ActorDensityLimit: model = $(par.model) is unknown")
    end

    return actor
end

function _finalize(actor::ActorDensityLimit)
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
    density_greenwald(dd::IMAS.dd, par::FUSEparameters__ActorDensityLimit, lim::IMAS.stability__limit)

Standard limit in density using IMAS greenwald fraction
Model Formulation: f_{GW,IMAS} < 1.0
Citation: 
"""
function density_greenwald(dd::IMAS.dd, par::FUSEparameters__ActorDensityLimit, lim::IMAS.stability__limit)
    lim.model.name = "Standard::IMAS_Greenwald"
    lim.model.formula = "IMAS.greenwald_fraction < 1.0"

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    model_value = IMAS.greenwald_fraction(eqt, cp1d)
    target_value = 1.0

    lim.model.fraction = model_value / target_value
end