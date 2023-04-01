#= ============== =#
#  ActorBetaLimit  #
#= ============== =#

Base.@kwdef mutable struct FUSEparameters__ActorBetaLimit{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:Li, :Standard, :None], "-", "Model for the limit calculation"; default=:None)
    submodel::Switch{Symbol} = Switch(Symbol, [:A, :B, :C, :None], "-", "Submodel for the limit calculation"; default=:None)
end

mutable struct ActorBetaLimit <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorBetaLimit
    limit::IMAS.stability__limit
end

"""
    ActorBetaLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)


"""
function ActorBetaLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorBetaLimit(kw...)
    actor = ActorBetaLimit(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorBetaLimit(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit; kw...)
    logging_actor_init(ActorBetaLimit)
    par = par(kw...)
    limit = resize!(dd.stability.limit, "name" => "Beta limit")
    ActorBetaLimit(dd, par, limit)
end

"""
    step(actor::ActorBetaLimit)

Runs ActorBetaLimit to evaluate the beta limit for the given equilibrium
"""
function _step(actor::ActorBetaLimit)
    dd = actor.dd
    par = actor.par
    limit = actor.limit

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal

    if par.model == :None 
        logging(Logging.Error, :actors, "ActorBetaLimit: limit check disabled")

    elseif par.model == :Standard
        model_value = beta_normal # This is the actual model being applied
        if par.submodel == :None # Fail
            error("ActorBetaLimit: model = $(par.model):$(par.submodel) is not implemented")
        elseif par.submodel == :A #Troyon Limit
            target_value = beta_standard_1(dd, par, limit)
        elseif par.submodel == :B #Classical Limit
            target_value = beta_standard_2(dd, par, limit)
        else
            error("ActorBetaLimit: model = $(par.model):$(par.submodel) is unknown")
        end

    elseif par.model == :Li
        #actor.val = model2(eqt)
    else
        error("ActorBetaLimit: model = $(par.model) is unknown")
    end

    limit.model.fraction = model_value / target_value
    return actor
end

function _finalize(actor::ActorBetaLimit)
    limit = actor.limit
    if limit.model.fraction  < 1
        limit.cleared = 1
    else
        limit.cleared = 0 
    end

    return actor
end


##### Limit Model Functions #####

# function model_standard(dd, par, limit)
#     dd = actor.dd
#     par = actor.par
#     limit = actor.limit

#     beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
#     model_value = beta_normal

# end

function beta_standard_1(dd, par, limit)
    limit.model.name = "Standard::Troyon"
    limit.model.formula = "Beta_{N} < 2.8"
    target_value = limit.model.limit = 2.8
    return target_value
end

function beta_standard_2(dd, par, limit)
    limit.model.name = "Standard::Classical"
    limit.model.formula = "Beta_{N} < 3.5"
    target_value = limit.model.limit = 3.5
    return target_value
end

# function model_name_number()
#     limit.model.name = "NAME"
#     #limit.model.description = "LONG DESCRIPTION"
#     limit.model.formula = "GENERAL FORMULA"
#     target_value = 1.0 # value 
#     actual_value = formula_result
#     fraction = actual_value / target_value
#     return fraction
# end

