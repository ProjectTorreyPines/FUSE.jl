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
    lim::IMAS.stability__model
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
    #lim = resize!(dd.stability.limit, "name" => "Density limit")
    model_index = 999
    lim = resize!(dd.stability.model, "identifier.index" => model_index)
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

    if par.model == :None 
        logging(Logging.Error, :actors, "ActorDensityLimit: limit check disabled")
    elseif par.model == :Standard
        density_standard_a(dd, par, lim)
    else
        error("ActorDensityLimit: model = $(par.model) is unknown")
    end

    return actor
end

function _finalize(actor::ActorDensityLimit)
    # lim = actor.lim
    # if lim.fraction  < 1
    #     lim.cleared = 1
    # else
    #     lim.cleared = 0 
    # end
    return actor
end

##### Limit Model Functions #####

"""
    density_standard(dd::IMAS.dd)

Standard limit in density
Model Formulation: f_GW < 1
Citation: 
"""
function density_standard(dd::IMAS.dd)
 
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    model_value = IMAS.greenwald_fraction(eqt, cp1d)

    return model_value
end

"""
    density_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorDensityLimit, lim::IMAS.stability__model)

Standard limit in density using IMAS greenwald fraction
Model Formulation: f_{GW,IMAS} < 1.0
Citation: 
"""
function density_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorDensityLimit, lim::IMAS.stability__model)
    lim.identifier.name = "Standard::IMAS_Greenwald"
    lim.identifier.description = "IMAS.greenwald_fraction < 1.0"

    model_value = density_standard(dd)
    target_value = 1.0

    lim.fraction = model_value / target_value
end