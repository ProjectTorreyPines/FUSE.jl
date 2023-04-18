##### SPECIAL CASES #####

function force_fail(dd::IMAS.dd, model::IMAS.stability__model)
    error(raw"¯\_(ツ)_/¯") #I'll change this eventually
end

function run_stability_models(dd::IMAS.dd, model_name::Symbol)
    return run_stability_models(dd, [model_name])
end

function run_stability_models(dd::IMAS.dd, model_names::Vector{Symbol})
    for model_name in model_names
        eval(model_name)(dd)
    end
end

##### MODEL COLLECTIONS #####
function default(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:default]
    model = resize!(dd.stability.collection, "identifier.index" => model_index)
    model.identifier.name = "Default"
    model.identifier.description = "Uses the default set of models"

    models = [:beta_troyon_1984, :model_201, :model_301, :model_401]
    run_stability_models(dd, models)
end

function beta_limits(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:beta_limits]
    model = resize!(dd.stability.collection, "identifier.index" => model_index)
    model.identifier.name = "Beta Limits"
    model.identifier.description = "Checks all beta limit models"

    models = [:beta_troyon_1984, :beta_troyon_1985, :beta_tuda_1985, :beta_bernard_1983, :beta_model_105]
    run_stability_models(dd, models)
end

function current_limits(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:current_limits]
    model = resize!(dd.stability.collection, "identifier.index" => model_index)
    model.identifier.name = "Current Limits"
    model.identifier.description = "Checks all current limit models"

    models = [:model_201]
    run_stability_models(dd, models)
end

function density_limits(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:density_limits]
    model = resize!(dd.stability.collection, "identifier.index" => model_index)
    model.identifier.name = "Density Limits"
    model.identifier.description = "Checks all density limit models"

    models = [:model_301]
    run_stability_models(dd, models)
end


##### BETA LIMIT MODELS #####

"""
    beta_troyon_1984(dd::IMAS.dd)

Limit in normalized beta using Troyon scaling

Model Formulation: Beta_{N} < 3.5

Citation:  F Troyon et al 1984 Plasma Phys. Control. Fusion 26 209
"""
function beta_troyon_1984(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:beta_troyon_1984]
    model = resize!(dd.stability.model, "identifier.index" => model_index)
    model.identifier.name = "Troyon 1984"
    model.identifier.description = "Beta_{N} < 3.5"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 3.5

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_troyon_1985(dd::IMAS.dd)

Limit in normalized beta using classical scaling using combined kink and ballooning stability

Model Formulation: Beta_{N} < 2.8

Citation: 
"""
function beta_troyon_1985(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:beta_troyon_1985]
    model = resize!(dd.stability.model, "identifier.index" => model_index)
    model.identifier.name = "Troyon 1985"
    model.identifier.description = "Beta_{N} < 2.8"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 2.8

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_tuda_1985(dd::IMAS.dd)

Limit in beta_normal using classical scaling using only kink stability

Model Formulation: Beta_{N} < 3.2

Citation: 
"""
function beta_tuda_1985(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:beta_tuda_1985]
    model = resize!(dd.stability.model, "identifier.index" => model_index)
    model.identifier.name = "Tuda 1985"
    model.identifier.description = "Beta_{N} < 3.2"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 3.2

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_bernard1983(dd::IMAS.dd)

Limit in normalized beta using classical scaling using only ballooning stability

Model Formulation: Beta_{N} < 2.8

Citation: 
"""
function beta_bernard_1983(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:beta_bernard_1983]
    model = resize!(dd.stability.model, "identifier.index" => model_index)
    model.identifier.name = "Bernard 1983"
    model.identifier.description = "Beta_{N} < 4.4"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 4.4

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_betali_a(dd::IMAS.dd)

Modern limit in normlaized beta normalized by plasma inductance

Model Formulation: beta_{N} / li < C_{beta}

Citation: 
"""
function model_105(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:model_105]
    model = resize!(dd.stability.model, "identifier.index" => model_index)
    model.identifier.name = "BetaLi::Troyon"
    model.identifier.description = "Beta_{N} / Li < 4.0"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    plasma_inductance = dd.equilibrium.time_slice[].global_quantities.li_3

    model_value = beta_normal / plasma_inductance
    target_value = 4.4

    @ddtime(model.fraction = model_value / target_value)
end


##### CURRENT LIMIT MODELS #####

"""
    model_201(dd::IMAS.dd)

Standard limit in edge current via the safety factor

Model Formulation: q95 < 2

Citation: 
"""
function model_201(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:model_201]
    model = resize!(dd.stability.model, "identifier.index" => model_index)
    model.identifier.name = "Standard::q95"
    model.identifier.description = "q_95 > 2.0"

    q95 = dd.equilibrium.time_slice[].global_quantities.q_95
    model_value = 1 / abs(q95)
    target_value = 0.5

    @ddtime(model.fraction = model_value / target_value)
end


##### DENSITY LIMIT MODELS #####

"""
    model_301(dd::IMAS.dd)

Standard limit in density using IMAS greenwald fraction

Model Formulation: f_{GW,IMAS} < 1.0

Citation: 
"""
function model_301(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:model_301]
    model = resize!(dd.stability.model, "identifier.index" => model_index)
    model.identifier.name = "Standard::IMAS_Greenwald"
    model.identifier.description = "IMAS.greenwald_fraction < 1.0"

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    model_value = IMAS.greenwald_fraction(eqt, cp1d)
    target_value = 1.0

    @ddtime(model.fraction = model_value / target_value)
end


##### SHAPE LIMIT MODELS #####

function model_401(dd::IMAS.dd)
    model_index = IMAS.name_2_index(dd.stability.model)[:model_401]
    model = resize!(dd.stability.model, "identifier.index" => model_index)
    model.identifier.name = "Standard elongation limit"
    model.identifier.description = "elongation < IMAS.elongation_limit"

    eqt = dd.equilibrium.time_slice[]
    model_value = dd.equilibrium.time_slice[].boundary.elongation
    target_value = IMAS.elongation_limit(eqt)

    @ddtime(model.fraction = model_value / target_value)

end

##### QUESTIONS #####
# 2) Add actual links to references in the docstring?
#    > How to do this? 
#    > Can we clone the docstring to the id.description?
# 4) Best way to handle inputs to keep compatability with function Dict?
#    > new mutable struct? function params? 
