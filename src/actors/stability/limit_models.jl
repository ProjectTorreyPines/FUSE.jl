##### SPECIAL CASES #####

function force_fail(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model) 
    error(raw"¯\_(ツ)_/¯") #I'll change this eventually
end

##### MODEL COLLECTIONS #####

function default(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__collection)
    model.identifier.name = "Default"
    model.identifier.description = "Uses the default set of models"

    logging(Logging.Error, :actors, "ActorStabilityLimits: default model may not be sufficent check on stability.")

    par.model_ids = [:troyon_1984, :model_201, :model_301, :model_401]
    step(ActorStabilityLimits(dd, par))
end

function beta_limits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__collection)
    model.identifier.name = "Beta Limits"
    model.identifier.description = "Checks all beta limit models"

    par.model_ids = [:troyon_1984, :troyon_1985, :tuda_1985, :bernard_1983, :model_105]
    step(ActorStabilityLimits(dd, par))
end

function current_limits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__collection)
    model.identifier.name = "Current Limits"
    model.identifier.description = "Checks all current limit models"

    par.model_ids = [:model_201]
    step(ActorStabilityLimits(dd, par))
end

function density_limits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__collection)
    model.identifier.name = "Density Limits"
    model.identifier.description = "Checks all density limit models"

    par.model_ids = [:model_301]
    step(ActorStabilityLimits(dd, par))
end


##### BETA LIMIT MODELS #####

"""
    beta_troyon_1984(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Limit in normalized beta using Troyon scaling
Model Formulation: Beta_{N} < 3.5
Citation:  F Troyon et al 1984 Plasma Phys. Control. Fusion 26 209
"""
#funciton model_101(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
function beta_troyon_1984(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Troyon 1984"
    model.identifier.description = "Beta_{N} < 3.5"
 
    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal 
    target_value = 3.5

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_troyon_1985(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Limit in normalized beta using classical scaling using combined kink and ballooning stability
Model Formulation: Beta_{N} < 2.8
Citation: 
"""
function beta_troyon_1985(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Troyon 1985"
    model.identifier.description = "Beta_{N} < 2.8"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal 
    target_value = 2.8

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_tuda_1985(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Limit in beta_normal using classical scaling using only kink stability
Model Formulation: Beta_{N} < 3.2
Citation: 
"""
function beta_tuda_1985(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Tuda 1985"
    model.identifier.description = "Beta_{N} < 3.2"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal 
    target_value = 3.2

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_bernard1983(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Limit in normalized beta using classical scaling using only ballooning stability
Model Formulation: Beta_{N} < 2.8
Citation: 
"""
function beta_bernard_1983(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Bernard 1983"
    model.identifier.description = "Beta_{N} < 4.4"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal 
    target_value = 4.4

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_betali_a(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Modern limit in normlaized beta normalized by plasma inductance
Model Formulation: beta_{N} / li < C_{beta}
Citation: 
"""
function model_105(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "BetaLi::Troyon"
    model.identifier.description = "Beta_{N} / Li < 4.0"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    plasma_inductance =  dd.equilibrium.time_slice[].global_quantities.li_3
    
    model_value = beta_normal / plasma_inductance
    target_value = 4.4

    @ddtime(model.fraction = model_value / target_value)
end


##### CURRENT LIMIT MODELS #####

"""
    model_201(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Standard limit in edge current via the safety factor
Model Formulation: q95 < 2
Citation: 
"""
function model_201(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Standard::q95"
    model.identifier.description = "q_95 > 2.0"

    q95 = dd.equilibrium.time_slice[].global_quantities.q_95 
    model_value = 1/abs(q95)
    target_value = 0.5

    @ddtime(model.fraction = model_value / target_value)
end


##### DENSITY LIMIT MODELS #####

"""
    model_301(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Standard limit in density using IMAS greenwald fraction
Model Formulation: f_{GW,IMAS} < 1.0
Citation: 
"""
function model_301(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Standard::IMAS_Greenwald"
    model.identifier.description = "IMAS.greenwald_fraction < 1.0"

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    model_value = IMAS.greenwald_fraction(eqt, cp1d)
    target_value = 1.0

    @ddtime(model.fraction = model_value / target_value)
end


##### SHAPE LIMIT MODELS #####

function model_401(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Standard elongation limit"
    model.identifier.description = "elongation < IMAS.elongation_limit"

    eqt = dd.equilibrium.time_slice[]
    model_value = dd.equilibrium.time_slice[].boundary.elongation
    target_value = IMAS.elongation_limit(eqt)

    @ddtime(model.fraction = model_value / target_value)

end

##### MAPPING DICTIONARY #####
# Is there a better way to do this? 

const limit_models = Dict(
    0 => force_fail,
    1 => default,
    11 => beta_limits,
    12 => current_limits,
    13 => density_limits,
    101 => beta_troyon_1984,
    102 => beta_troyon_1985,
    103 => beta_tuda_1985,
    104 => beta_bernard_1983,
    105 => model_105,
    201 => model_201,
    301 => model_301,
    401 => model_401
)


##### QUESTIONS #####
# 1) Putting the models in a constant Dictionary of functions
#    > should match the constant dataset
#    > Is there a better way to link them? 
#    > something more automatic?
# 2) Add actual links to references in the docstring?
#    > How to do this? 
#    > Can we clone the docstring to the id.description?
# 3) Best way to implement collections?
#    > Just run all the model functions?
#    >>> Changed `dd` to make this work better
#    > Run new actors for each? 
# 4) Best way to handle inputs to keep compatability with function Dict?
#    > new mutable struct? function params? 



 









