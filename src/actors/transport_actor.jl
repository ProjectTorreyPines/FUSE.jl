using NumericalIntegration
import AD_TAUENN

function init(cs::IMAS.core_sources, eq::IMAS.equilibrium; Paux_e::Real, Paux_i::Real, ngrid::Int = 51)
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

function init(
    cp::IMAS.core_profiles,
    eq::IMAS.equilibrium,
    summary::IMAS.summary;
    ne_ped::Real,
    ne_peaking::Real,
    Te_ped::Real,
    Te_peaking::Real,
    w_ped::Real,
    zeff::Real,
    Paux::Real,
    T_ratio::Real = 1.0,
    n_points::Int = 101,
)
    cpt = resize!(cp.profiles_1d)

    cpt.grid.rho_tor_norm = LinRange(0, 1, n_points)
    cpt.zeff = ones(n_points) .* zeff
    cpt.rotation_frequency_tor_sonic = 5e3 * abs(Paux / 1e6 * 1.0 + 0.5) .* (1.0 .- cpt.grid.rho_tor_norm)

    # Set ions
    resize!(cpt.ion, 2)
    cpt.ion[1].label = "D"
    resize!(cpt.ion[1].element, 1)
    cpt.ion[1].element[1].z_n = 1
    cpt.ion[1].element[1].a = 2
    resize!(cpt.ion[2].element, 1)
    cpt.ion[2].label = "C"
    cpt.ion[2].element[1].z_n = 6
    cpt.ion[2].element[1].a = 12

    # pedestal
    @ddtime summary.local.pedestal.n_e.value = ne_ped
    @ddtime summary.local.pedestal.position.rho_tor_norm = 1 - w_ped
    @ddtime summary.local.pedestal.zeff.value = zeff

    # Set densities
    ne_core = ne_peaking * ne_ped
    cpt.electrons.density = AD_TAUENN.Hmode_profiles(0.5 * ne_ped, ne_ped, ne_core, n_points, 0.9, 0.9, w_ped)

    zimp1 = 6.0
    niFraction = zeros(2)
    niFraction[2] = (zeff - 1.0) / (zimp1 * (zimp1 - 1.0))
    niFraction[1] = 1.0 - zimp1 * niFraction[2]
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

    # to be done as an expression
    # prof1d.pressure_thermal = 1.6e-19 .* cpt.electrons.density .* cpt.electrons.temperature
    # for i in 1:length(prof1d.ion)
    #     prof1d.pressure_thermal .+= 1.6e-19 .* prof1d.ion[i].density .* prof1d.ion[i].temperature
    # end
    return cp
end

#= ================ =#
#     TAUENN actor   #
#= ================ =#

mutable struct TaueNNactor <: AbstractActor
    dd::IMAS.dd
    rho_fluxmatch::Real
    eped_factor::Real
    temp_shape::Real
    temp_pedestal_ratio::Real
end

function TaueNNactor(dd::IMAS.dd; rho_fluxmatch = 0.5, eped_factor = 1.0, temp_shape = 1.8, temp_pedestal_ratio = 1.0)
    return TaueNNactor(dd, rho_fluxmatch, eped_factor, temp_shape, temp_pedestal_ratio)
end

function step(actor::TaueNNactor)
    return AD_TAUENN.tau_enn(actor.dd, actor.rho_fluxmatch, actor.eped_factor, actor.temp_shape, actor.temp_pedestal_ratio)
end