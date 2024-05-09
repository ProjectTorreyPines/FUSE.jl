function supported_stability_models()
    return "[:$(join(values(IMAS.index_2_name__stability__collection),", :")), :$(join(values(IMAS.index_2_name__stability__model),", :"))]"
end

function run_stability_models(dd::IMAS.dd, model_name::Symbol)
    return run_stability_models(dd, [model_name])
end

function run_stability_models(dd::IMAS.dd, model_names::Vector{Symbol})
    for model_name in model_names
        func = try
            eval(model_name)
        catch e
            if typeof(e) <: UndefVarError
                error("Unknown model `$(repr(model_name))`. Supported models are $(supported_stability_models())")
            else
                rethrow(e)
            end
        end
        func(dd)
    end
end

##### MODEL COLLECTIONS #####

function default_limits(dd::IMAS.dd)
    collection = resize!(dd.stability.collection, :default_limits)
    collection.identifier.name = "Default Limits"
    collection.identifier.description = "Uses the default set of models"

    models = [:beta_troyon_1984, :q95_gt_2, :gw_density, :κ_controllability]
    return run_stability_models(dd, models)
end

function beta_limits(dd::IMAS.dd)
    collection = resize!(dd.stability.collection, :beta_limits)
    collection.identifier.name = "Beta Limits"
    collection.identifier.description = "Checks all beta limit models"

    models = [:beta_troyon_1984, :beta_troyon_1985, :beta_tuda_1985, :beta_bernard_1983]
    return run_stability_models(dd, models)
end

function current_limits(dd::IMAS.dd)
    collection = resize!(dd.stability.collection, :current_limits)
    collection.identifier.name = "Current Limits"
    collection.identifier.description = "Checks all current limit models"

    models = [:q95_gt_2]
    return run_stability_models(dd, models)
end

function density_limits(dd::IMAS.dd)
    collection = resize!(dd.stability.collection, :density_limits)
    collection.identifier.name = "Density Limits"
    collection.identifier.description = "Checks all density limit models"

    models = [:gw_density]
    return run_stability_models(dd, models)
end


##### BETA LIMIT MODELS #####
"""
    beta_troyon_1984(dd::IMAS.dd)

Limit in normalized beta using Troyon scaling

Model formulation: `βn < 3.5`

Citation:  F Troyon et al 1984 Plasma Phys. Control. Fusion 26 209
"""
function beta_troyon_1984(dd::IMAS.dd)
    model = resize!(dd.stability.model, :beta_troyon_1984)
    model.identifier.name = "Troyon 1984"
    model.identifier.description = "βn < 3.5"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 3.5

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_troyon_1985(dd::IMAS.dd)

Limit in normalized beta using classical scaling using combined kink and ballooning stability

Model formulation: `βn < 2.8`

Citation:
"""
function beta_troyon_1985(dd::IMAS.dd)
    model = resize!(dd.stability.model, :beta_troyon_1985)
    model.identifier.name = "Troyon 1985"
    model.identifier.description = "βn < 2.8"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 2.8

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_tuda_1985(dd::IMAS.dd)

Limit in beta_normal using classical scaling using only kink stability

Model formulation: `βn < 3.2`

Citation:
"""
function beta_tuda_1985(dd::IMAS.dd)
    model = resize!(dd.stability.model, :beta_tuda_1985)
    model.identifier.name = "Tuda 1985"
    model.identifier.description = "βn < 3.2"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 3.2

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_bernard1983(dd::IMAS.dd)

Limit in normalized beta using classical scaling using only ballooning stability

Model formulation: `βn < 2.8`

Citation:
"""
function beta_bernard_1983(dd::IMAS.dd)
    model = resize!(dd.stability.model, :beta_bernard_1983)
    model.identifier.name = "Bernard 1983"
    model.identifier.description = "βn < 4.4"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 4.4

    @ddtime(model.fraction = model_value / target_value)
end

"""
    beta_betali_a(dd::IMAS.dd)

Modern limit in normlaized beta normalized by plasma inductance

Model formulation: `βn / li < C_{beta}`

Citation:
"""
function beta_model_105(dd::IMAS.dd)
    model = resize!(dd.stability.model, :beta_model_105)
    model.identifier.name = "BetaLi::Troyon"
    model.identifier.description = "βn / Li < 4.0"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    plasma_inductance = dd.equilibrium.time_slice[].global_quantities.li_3

    model_value = beta_normal / plasma_inductance
    target_value = 4.4

    @ddtime(model.fraction = model_value / target_value)
end


##### CURRENT LIMIT MODELS #####

"""
    q95_gt_2(dd::IMAS.dd)

Standard limit in edge current via the safety factor

Model formulation: `q95 > 2.0`

Citation:
"""
function q95_gt_2(dd::IMAS.dd)
    model = resize!(dd.stability.model, :q95_gt_2)
    model.identifier.name = "Standard::q95"
    model.identifier.description = "q_95 > 2.0"

    q95 = dd.equilibrium.time_slice[].global_quantities.q_95
    model_value = abs(q95)
    target_value = 2.0

    @ddtime(model.fraction = target_value / model_value)
end

"""
    q08_gt_2(dd::IMAS.dd)

Limit in edge current via the safety factor `q(rho=0.8) > 2.0`

Model formulation: `q(rho=0.8) > 2.0`
"""
function q08_gt_2(dd::IMAS.dd)
    model = resize!(dd.stability.model, :q08_gt_2)
    model.identifier.name = "q(rho=0.8) > 2.0"
    model.identifier.description = "q(rho=0.8) > 2.0"

    rho_eq = dd.equilibrium.time_slice[].profiles_1d.rho_tor_norm
    q_08 = abs.(dd.equilibrium.time_slice[].profiles_1d.q)[argmin(abs.(rho_eq .- 0.8))]

    model_value = abs(q_08)
    target_value = 2.0

    @ddtime(model.fraction = target_value / model_value)
end

##### DENSITY LIMIT MODELS #####

"""
    gw_density(dd::IMAS.dd)

Standard limit in density using IMAS greenwald fraction

Model formulation: `f_{GW,IMAS} < 1.0`

Citation:
"""
function gw_density(dd::IMAS.dd)
    model = resize!(dd.stability.model, :gw_density)
    model.identifier.name = "Standard::IMAS_Greenwald"
    model.identifier.description = "IMAS.greenwald_fraction < 1.0"

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    model_value = IMAS.greenwald_fraction(eqt, cp1d)
    target_value = 1.0

    @ddtime(model.fraction = model_value / target_value)
end


##### SHAPE LIMIT MODELS #####

function κ_controllability(dd::IMAS.dd)
    model = resize!(dd.stability.model, :κ_controllability)
    model.identifier.name = "Standard elongation limit"
    model.identifier.description = "elongation < IMAS.elongation_limit"

    eqt = dd.equilibrium.time_slice[]
    model_value = dd.equilibrium.time_slice[].boundary.elongation
    target_value = IMAS.elongation_limit(eqt)

    @ddtime(model.fraction = model_value / target_value)
end