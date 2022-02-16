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


