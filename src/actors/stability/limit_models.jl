const supported_limit_models = Symbol[]
const default_limit_models = Symbol[]

##### VERTICAL STABILITY #####

function vertical_stability(dd::IMAS.dd, act::ParametersAllActors)
    ActorVerticalStability(dd, act)
    mhd = dd.mhd_linear.time_slice[]

    # Vertical stability margin > 0.15 for stability
    model = resize!(dd.limits.model, "identifier.name" => "Vertical Stability margin")
    model.identifier.description = "Vertical stability margin > 0.15 for stability"
    index = findfirst(mode -> mode.perturbation_type.name == "m_s" && mode.n_tor == 0, mhd.toroidal_mode)
    if index !== nothing
        mode = mhd.toroidal_mode[index]
        model_value = mode.stability_metric
        target_value = 0.15
        @ddtime(model.fraction = target_value / model_value)
    else
        @ddtime(model.fraction = 0.0)
    end

    # Normalized vertical growth rate, < 10 for stability
    model = resize!(dd.limits.model, "identifier.name" => "Vertical Stability growth rate")
    model.identifier.description = "Normalized vertical growth rate < 10 for stability"
    index = findfirst(mode -> mode.perturbation_type.name == "γτ" && mode.n_tor == 0, mhd.toroidal_mode)
    if index !== nothing
        mode = mhd.toroidal_mode[index]
        model_value = mode.stability_metric
        target_value = 10.0
        @ddtime(model.fraction = model_value / target_value)
    else
        @ddtime(model.fraction = 0.0)
    end
end
push!(supported_limit_models, :vertical_stability)
push!(default_limit_models, :vertical_stability)

##### BETA LIMIT MODELS #####

function beta_troyon_nn(dd::IMAS.dd, act::ParametersAllActors)
    ActorTroyonBetaNN(dd, act)
    mhd = dd.mhd_linear.time_slice[]
    for n in 1:3
        model = resize!(dd.limits.model, "identifier.name" => "BetaTroyonNN n=$(n)")
        model.identifier.description = "βn < BetaTroyonNN n=$(n)"

        beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
        model_value = beta_normal

        index = findfirst(mode -> mode.perturbation_type.name == "Troyon no-wall" && mode.n_tor == n, mhd.toroidal_mode)
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
    target_value = 4.0

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

"""
    edge_collisionality(dd::IMAS.dd)

L-mode density limit based on edge collisionality 

Model formulation: `ν_limit_edge = 3.0 * β_T_edge ^ -0.41`

Sources: https://arxiv.org/pdf/2406.18442 (Maris et al, 2024) and https://iopscience.iop.org/article/10.1088/1741-4326/abdb91 (Verdoolaege et al, 2024)
"""
function edge_collisionality(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    model = resize!(dd.stability.model, :edge_collisionality)
    model.identifier.name = "Edge collisionality limit"
    model.identifier.description = "edge ν* < 3.0 * edge βT^-0.41"

    R0 = eqt.global_quantities.vacuum_toroidal_field.r0
    B0 = abs(eqt.global_quantities.vacuum_toroidal_field.b0)

    Rgeo = (eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]) / 2.0
    Btvac = B0 * R0 / Rgeo

    # edge is defined as all points between rho = 0.85 and rho = 0.95
    edge = findall(x -> 0.85 <= x <= 0.95, dd.core_profiles.profiles_1d[].grid.rho_tor_norm)
    edge_density = sum(cp1d.electrons.density[edge]) / length(edge) # 1/m^3
    edge_temperature = sum(cp1d.electrons.temperature[edge]) / length(edge) # eV 
    
    # Equation 4, with the addition of a factor of e to convert temperature from eV into Joules 
    beta_tor_edge = 2.0 * edge_density * constants.e * edge_temperature * 1e2 / (Btvac^2 / (2.0 * constants.μ_0))

    # Definition of loglam from Verdoolaege - Updated ITPA global confinement database, page 10
    loglam = 30.9 .- log.(cp1d.electrons.density.^(1/2) ./ cp1d.electrons.temperature)
    loglam_edge = sum(loglam[edge]) / length(edge)
    epsilon = eqt.boundary.minor_radius / R0 
    kappa = eqt.boundary.elongation 
    ip = eqt.global_quantities.ip

    # Definition of nu_star from Maris Equation 2, page 7 
    first_term = constants.e^2 * loglam_edge / (3^(3/2)* 2 * pi * constants.ϵ_0^2)
    second_term = edge_density / (edge_temperature^2)
    qcyl = (2*pi / constants.μ_0) * ((Btvac * kappa * eqt.boundary.minor_radius^2) / (ip * R0)) 
    third_term = (qcyl * R0) / (epsilon^(3/2))

    ν_star_edge = first_term * second_term * third_term

    model_value = ν_star_edge
    target_value = 3.0 * beta_tor_edge^-0.41

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