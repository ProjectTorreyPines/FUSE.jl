function run_limits_models(dd::IMAS.dd, act::ParametersAllActors, model_name::Symbol)
    return run_limits_models(dd, act, [model_name])
end

function run_limits_models(dd::IMAS.dd, act::ParametersAllActors, model_names::Vector{Symbol})
    for model_name in model_names
        func = try
            eval(model_name)
        catch e
            if typeof(e) <: UndefVarError
                error("Unknown model `$(repr(model_name))`. Supported models are $(supported_limit_models)")
            else
                rethrow(e)
            end
        end
        func(dd, act)
    end
end

supported_limit_models = Symbol[]
default_limit_models = Symbol[]

##### BETA LIMIT MODELS #####

function beta_troyon_nn(dd::IMAS.dd, act::ParametersAllActors)
    ActorTroyonBetaNN(dd, act)
    mhd = dd.mhd_linear.time_slice[]

    for n in 1:3
        model = resize!(dd.limits.model, "identifier.name" => "BetaTroyonNN__n$(n)")
        model.identifier.description = "βn < BetaTroyonNN n=$(n)"

        beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
        model_value = beta_normal

        index = findfirst(mode -> mode.perturbation_type.name == "Troyon no-wall" && mode.n_tor==n, mhd.toroidal_mode)
        target_value = mhd.toroidal_mode[index].stability_metric

        @ddtime(model.fraction = model_value / target_value)
    end
end
push!(supported_limit_models, :beta_troyon_nn)
push!(default_limit_models, :beta_troyon_nn)

"""
    beta_troyon_1984(dd::IMAS.dd, act::ParametersAllActors)

Limit in normalized beta using Troyon scaling

Model formulation: `βn < 3.5`

Citation:  F Troyon et al 1984 Plasma Phys. Control. Fusion 26 209
"""
function beta_troyon_1984(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "Troyon 1984")
    model.identifier.description = "βn < 3.5"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 3.5

    @ddtime(model.fraction = model_value / target_value)
end
push!(supported_limit_models, :beta_troyon_1984)

"""
    beta_troyon_1985(dd::IMAS.dd, act::ParametersAllActors)

Limit in normalized beta using classical scaling using combined kink and ballooning stability

Model formulation: `βn < 2.8`

Citation:
"""
function beta_troyon_1985(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "Troyon 1985")
    model.identifier.description = "βn < 2.8"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 2.8

    @ddtime(model.fraction = model_value / target_value)
end
push!(supported_limit_models, :beta_troyon_1985)

"""
    beta_tuda_1985(dd::IMAS.dd, act::ParametersAllActors)

Limit in beta_normal using classical scaling using only kink stability

Model formulation: `βn < 3.2`

Citation:
"""
function beta_tuda_1985(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "Tuda 1985")
    model.identifier.description = "βn < 3.2"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 3.2

    @ddtime(model.fraction = model_value / target_value)
end
push!(supported_limit_models, :beta_tuda_1985)

"""
    beta_bernard1983(dd::IMAS.dd, act::ParametersAllActors)

Limit in normalized beta using classical scaling using only ballooning stability

Model formulation: `βn < 2.8`

Citation:
"""
function beta_bernard_1983(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "Bernard 1983")
    model.identifier.description = "βn < 2.8"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal
    target_value = 2.8

    @ddtime(model.fraction = model_value / target_value)
end
push!(supported_limit_models, :beta_bernard_1983)

"""
    beta_betali_a(dd::IMAS.dd, act::ParametersAllActors)

Modern limit in normlaized beta normalized by plasma inductance

Model formulation: `βn / li < C_{beta}`

Citation:
"""
function beta_betali_a(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "BetaLi::Troyon")
    model.identifier.description = "βn / Li < 4.4"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    plasma_inductance = dd.equilibrium.time_slice[].global_quantities.li_3

    model_value = beta_normal / plasma_inductance
    target_value = 4.4

    @ddtime(model.fraction = model_value / target_value)
end
push!(supported_limit_models, :beta_betali_a)


##### CURRENT LIMIT MODELS #####

"""
    q95_gt_2(dd::IMAS.dd, act::ParametersAllActors)

Standard limit in edge current via the safety factor

Model formulation: `q95 > 2.0`

Citation:
"""
function q95_gt_2(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "q95 > 2.0")
    model.identifier.description = "q(rho=0.95) > 2.0"

    q95 = dd.equilibrium.time_slice[].global_quantities.q_95
    model_value = abs(q95)
    target_value = 2.0

    @ddtime(model.fraction = target_value / model_value)
end
push!(supported_limit_models, :q95_gt_2)
push!(default_limit_models, :q95_gt_2)

"""
    q80_gt_2(dd::IMAS.dd, act::ParametersAllActors)

Limit in edge current via the safety factor `q(rho=0.8) > 2.0`

Model formulation: `q(rho=0.8) > 2.0`
"""
function q80_gt_2(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "q80 > 2.0")
    model.identifier.description = "q(rho=0.8) > 2.0"

    rho_eq = dd.equilibrium.time_slice[].profiles_1d.rho_tor_norm
    q_08 = abs.(dd.equilibrium.time_slice[].profiles_1d.q)[argmin(abs.(rho_eq .- 0.8))]

    model_value = abs(q_08)
    target_value = 2.0

    @ddtime(model.fraction = target_value / model_value)
end
push!(supported_limit_models, :q80_gt_2)

##### DENSITY LIMIT MODELS #####

"""
    gw_density(dd::IMAS.dd, act::ParametersAllActors)

Standard limit in density using IMAS greenwald fraction

Model formulation: `f_{GW,IMAS} < 1.0`

Citation:
"""
function gw_density(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "greenwald_fraction < 1.0")
    model.identifier.description = "greenwald_fraction < 1.0"

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    model_value = IMAS.greenwald_fraction(eqt, cp1d)
    target_value = 1.0

    @ddtime(model.fraction = model_value / target_value)
end
push!(supported_limit_models, :gw_density)
push!(default_limit_models, :gw_density)

##### SHAPE LIMIT MODELS #####

function κ_controllability(dd::IMAS.dd, act::ParametersAllActors)
    model = resize!(dd.limits.model, "identifier.name" => "κ < elongation_limit")
    model.identifier.description = "elongation < IMAS.elongation_limit"

    eqt = dd.equilibrium.time_slice[]
    model_value = dd.equilibrium.time_slice[].boundary.elongation
    target_value = IMAS.elongation_limit(eqt)

    @ddtime(model.fraction = model_value / target_value)
end
push!(supported_limit_models, :κ_controllability)
push!(default_limit_models, :κ_controllability)