using NumericalIntegration
import AD_TAUENN

function init_core_sources(dd::IMAS.dd; kwargs...)
    init_core_sources(dd.core_sources, dd.equilibrium; kwargs...)
end

function init_core_sources(cs::IMAS.core_sources, eq::IMAS.equilibrium; Paux_e::Real, Paux_i::Real, ngrid::Int = 51)
    empty!(cs)

    resize!(cs.source, 1)
    cs.source[1].identifier.name = "arb"
    cs.source[1].identifier.index = 901
    cs.source[1].identifier.description = "Arbitrary source from FUSE transport initialization"

    resize!(cs.source[1].profiles_1d)
    cs1d = cs.source[1].profiles_1d[]
    cs1d.grid.rho_tor_norm = rho = LinRange(0.0, 1.0, ngrid)

    rho_eq = eq.time_slice[].profiles_1d.rho_tor_norm
    cs1d.grid.volume = IMAS.interp(rho_eq, eq.time_slice[].profiles_1d.volume)[rho]

    auxHeatingProfile = exp.(-4.0 * rho)
    pow_prof = cumul_integrate(cs1d.grid.volume, auxHeatingProfile)
    pow_prof = pow_prof ./ pow_prof[end]

    cs1d.electrons.power_inside = pow_prof .* Paux_e
    cs1d.total_ion_power_inside = pow_prof .* Paux_i

    return cs
end

function init_core_profiles(dd::IMAS.dd; kwargs...)
    init_core_profiles(dd.core_profiles, dd.equilibrium, dd.summary; kwargs...)
end

function init_core_profiles(
    cp::IMAS.core_profiles,
    eq::IMAS.equilibrium,
    summary::IMAS.summary;
    ne_ped::Real,
    ne_peaking::Real,
    Te_ped::Real,
    Te_peaking::Real,
    w_ped::Real,
    zeff::Real,
    bulk::Symbol,
    impurity::Symbol,
    Paux::Real,
    T_ratio::Real = 1.0,
    n_points::Int = 101
)
    cpt = resize!(cp.profiles_1d)

    cpt.grid.rho_tor_norm = LinRange(0, 1, n_points)
    cpt.zeff = ones(n_points) .* zeff
    cpt.rotation_frequency_tor_sonic = 5e3 * abs(Paux / 1e6 * 1.0 + 0.5) .* (1.0 .- cpt.grid.rho_tor_norm)

    # Set ions
    ion = resize!(cpt.ion, "label" => String(bulk))
    fill!(ion, IMAS.ion_element(bulk))
    @assert ion.element[1].z_n == 1 "Bulk ion must be a Hydrogen isotope [:H, :D, :DT, :T]"
    ion = resize!(cpt.ion, "label" => String(impurity))
    fill!(ion, IMAS.ion_element(impurity))

    # pedestal
    @ddtime summary.local.pedestal.n_e.value = ne_ped
    @ddtime summary.local.pedestal.position.rho_tor_norm = 1 - w_ped
    @ddtime summary.local.pedestal.zeff.value = zeff

    # Set densities
    ne_core = ne_peaking * ne_ped
    cpt.electrons.density = AD_TAUENN.Hmode_profiles(0.5 * ne_ped, ne_ped, ne_core, n_points, 0.9, 0.9, w_ped)

    zimp1 = IMAS.ion_element(impurity).element[1].z_n
    niFraction = zeros(2)
    niFraction[2] = (zeff - 1.0) / (zimp1 * (zimp1 - 1.0))
    niFraction[1] = 1.0 - zimp1 * niFraction[2]
    @assert all(niFraction .> 0.0) "zeff too high for the given bulk [$bulk] and impurity [$impurity] species"
    for i = 1:length(cpt.ion)
        cpt.ion[i].density = cpt.electrons.density .* niFraction[i]
    end

    # Set temperatures
    eqt = eq.time_slice[]
    betaN = eqt.global_quantities.beta_normal
    Bt = @ddtime eq.vacuum_toroidal_field.b0
    Ip = eqt.global_quantities.ip
    a = eqt.boundary.minor_radius

    Te_core = 10.0 * betaN * abs(Bt * (Ip / 1e6)) / a / (ne_core / 1e20) / (2.0 * 1.6e1 * 4.0 * pi * 1.0e-4)
    cpt.electrons.temperature = AD_TAUENN.Hmode_profiles(80.0, Te_ped, Te_core, n_points, Te_peaking, Te_peaking, w_ped)
    for i = 1:length(cpt.ion)
        cpt.ion[i].temperature = cpt.electrons.temperature ./ T_ratio
    end

    return cp
end

#= ================ =#
#     TAUENN actor   #
#= ================ =#

mutable struct TaueNNactor <: AbstractActor
    dd::IMAS.dd
    parameters::AD_TAUENN.TauennParameters
end

function TaueNNactor(dd::IMAS.dd; rho_fluxmatch = 0.6, eped_factor = 1.0, temp_shape = 1.8, temp_pedestal_ratio = 1.0, error=1E-2, use_tglfnn=true, kw...)
    parameters = AD_TAUENN.TauennParameters()
    parameters.eped_factor = eped_factor
    parameters.rho_fluxmatch = rho_fluxmatch
    parameters.temp_shape = temp_shape
    parameters.temp_pedestal_ratio = temp_pedestal_ratio
    parameters.use_tglfnn = use_tglfnn
    parameters.error = error
    for param in keys(kw)
        setfield!(parameters, param, kw[param])
    end
    return TaueNNactor(dd, parameters)
end

function step(actor::TaueNNactor; verbose = false)
    return AD_TAUENN.tau_enn(actor.dd, actor.parameters; verbose)
end