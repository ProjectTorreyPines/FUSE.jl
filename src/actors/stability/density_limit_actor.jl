#= ================= =#
#  ActorDensityLimit  #
#= ================= =#

Base.@kwdef mutable struct FUSEparameters__ActorDensityLimit{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:Standard, :None], "-", "Model for the limit calculation"; default=:None)
    submodel::Switch{Symbol} = Switch(Symbol, [:A, :B, :C, :None], "-", "Submodel for the limit calculation"; default=:None)
end

mutable struct ActorDensityLimit <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorDensityLimit
    limit::IMAS.stability__limit
end

"""
    ActorDensityLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)


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
    limit = resize!(dd.stability.limit, "name" => "Density limit")
    ActorDensityLimit(dd, par, limit)
end


"""
    step(actor::ActorDensityLimit)

Runs ActorDensityLimit to evaluate the Density limit for the given equilibrium
"""
function _step(actor::ActorDensityLimit)
    dd = actor.dd
    par = actor.par
    limit = actor.limit

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    if par.model == :None 
        logging(Logging.Error, :actors, "ActorDensityLimit: limit check disabled")
    elseif par.model == :Standard
        model_value = IMAS.greenwald_fraction(eqt, cp1d)
        target_value = density_standard_1(dd, par, limit)
    else
        error("ActorDensityLimit: model = $(par.model) is unknown")
    end

    limit.model.fraction = model_value / target_value
    return actor
end

function _finalize(actor::ActorDensityLimit)
    limit = actor.limit
    if limit.model.fraction  < 1
        limit.cleared = 1
    else
        limit.cleared = 0 
    end

    return actor
end

##### Limit Model Functions #####

function density_standard_1(dd, par, limit)
    limit.model.name = "Standard::IMAS_Greenwald"
    limit.model.formula = "IMAS.greenwald_fraction < 1.0"
    target_value = limit.model.limit = 1.0
    return target_value
end