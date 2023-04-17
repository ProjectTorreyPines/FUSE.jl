#= ============== =#
#  ActorBetaLimit  #
#= ============== =#

Base.@kwdef mutable struct FUSEparameters__ActorBetaLimit{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:BetaLi, :Standard, :None], "-", "Model for the limit calculation"; default=:None)
    submodel::Switch{Symbol} = Switch(Symbol, [:A, :B, :C, :D, :None], "-", "Submodel for the limit calculation"; default=:None)
end

mutable struct ActorBetaLimit <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorBetaLimit
    lim::IMAS.stability__model
end

"""
    ActorBetaLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates the emperical beta limit using various models
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

    @ddtime(dd.stability.time = dd.global_time)

    #model = par.model
    #model_index = IMAS.name_2_index(dd.stability.model)[model]
    model_index = 999
    lim = resize!(dd.stability.model, "identifier.index" => model_index)
    ActorBetaLimit(dd, par, lim)
end

"""
    step(actor::ActorBetaLimit)

Runs ActorBetaLimit to evaluate the beta limit for the given equilibrium
"""
function _step(actor::ActorBetaLimit)
    dd = actor.dd
    par = actor.par
    lim = actor.lim

    if par.model == :None 
        logging(Logging.Error, :actors, "ActorBetaLimit: limit check disabled")

    elseif par.model == :Standard
        if par.submodel == :None # Fail
            error("ActorBetaLimit: model = $(par.model):$(par.submodel) is not implemented")
        elseif par.submodel == :A #Troyon Limit
            beta_standard_a(dd, par, lim)
        elseif par.submodel == :B #Classical Limit
            beta_standard_b(dd, par, lim)
        elseif par.submodel == :C #Kink Only Limit
            beta_standard_c(dd, par, lim)
        elseif par.submodel == :D #Ballooning Only Limit
            beta_standard_d(dd, par, lim)
        else
            error("ActorBetaLimit: model = $(par.model):$(par.submodel) is unknown")
        end

    elseif par.model == :BetaLi
        if par.submodel == :None # Fail
            error("ActorBetaLimit: model = $(par.model):$(par.submodel) is not implemented")
        elseif par.submodel == :A #NAME Limit
            beta_betali_a(dd, par, lim)
        else
            error("ActorBetaLimit: model = $(par.model):$(par.submodel) is unknown")
        end
    else
        error("ActorBetaLimit: model = $(par.model) is unknown")
    end

    return actor
end

"""
    finalize(actor::ActorBetaLimit)

Writes ActorBetaLimit results to dd.stability
"""
function _finalize(actor::ActorBetaLimit)
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
    beta_standard(dd::IMAS.dd)

Standard limit in beta_normal
Model Formulation: beta_{N} < C_{beta}
Citation: 
"""
function beta_standard(dd::IMAS.dd)
    
    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    
    model_value = beta_normal 

    return model_value
end

"""
    beta_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)

Standard limit in beta_normal using Troyon scaling
Model Formulation: Beta_{N} < 3.5
Citation:  F Troyon et al 1984 Plasma Phys. Control. Fusion 26 209
"""
function beta_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)
    lim.identifier.name = "Standard::Troyon"
    lim.identifier.description = "Beta_{N} < 3.5"
    lim.identifier.index = 101
 
    model_value = beta_standard(dd)
    target_value = 3.5

    println(lim)

    @ddtime(lim.fraction = [model_value / target_value])
    #resize!(dd.stability.model, "fraction" => [fraction])
end

"""
    beta_standard_b(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)

Standard limit in beta_normal using classical scaling using combined kink and ballooning stability
Model Formulation: Beta_{N} < 2.8
Citation: 
"""
function beta_standard_b(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)
    lim.identifier.name = "Standard::Classical"
    lim.identifier.description = "Beta_{N} < 2.8"

    model_value = beta_standard(dd)
    target_value = 2.8

    @ddtime(lim.fraction = model_value / target_value)
end

"""
    beta_standard_c(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)

Standard limit in beta_normal using classical scaling using only kink stability
Model Formulation: Beta_{N} < 2.8
Citation: 
"""
function beta_standard_c(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)
    lim.identifier.name = "Standard::KinkOnly"
    lim.identifier.description = "Beta_{N} < 3.2"

    model_value = beta_standard(dd)
    target_value = 3.2

    @ddtime(lim.fraction = model_value / target_value)
end

"""
    beta_standard_d(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)

Standard limit in normalized beta using classical scaling using only ballooning stability
Model Formulation: Beta_{N} < 2.8
Citation: 
"""
function beta_standard_d(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)
    lim.identifier.name = "Standard::BallooningOnly"
    lim.identifier.description = "Beta_{N} < 4.4"

    model_value = beta_standard(dd)
    target_value = 4.4

    @ddtime(lim.fraction = model_value / target_value)
end

"""
    beta_betali(dd::IMAS.dd)

Modern limit in normlaized beta normalized by plasma inductance
Model Formulation: beta_{N} / li < C_{beta}
Citation: 
"""
function beta_betali(dd::IMAS.dd)
    
    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    plasma_inductance =  dd.equilibrium.time_slice[].global_quantities.li_3
    
    model_value = beta_normal / plasma_inductance

    return model_value
end

"""
    beta_betali_a(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)

Modern limit in normlaized beta normalized by plasma inductance
Model Formulation: beta_{N} / li < C_{beta}
Citation: 
"""
function beta_betali_a(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit, lim::IMAS.stability__model)
    lim.identifier.name = "BetaLi::Troyon"
    lim.identifier.description = "Beta_{N} / Li < 4.0"

    model_value = beta_betali(dd)
    target_value = 4.4

    @ddtime(lim.fraction = model_value / target_value)
end

# function model_name_number()
#     limit.lim.identifier.name = "NAME"
#     #limit.lim.description = "LONG DESCRIPTION"
#     limit.lim.identifier.description = "GENERAL FORMULA"
#     target_value = 1.0 # value 
#     actual_value = formula_result
#     fraction = actual_value / target_value
#     return fraction
# end

