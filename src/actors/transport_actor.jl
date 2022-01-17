


function init(cs::IMAS.core_sources, eq::IMAS.equilibrium, Paux_e::Real, Paux_i::Real)
    for time_index in 1:length(eq.time_slice)
        cs1d = cs.source[1].profiles_1d[time_index]
        cs1d = cs.source[1].index ...
        cs1d = cs.source[1].name ...
        cs1d.grid.rho_tor_norm = rho = eq.time_slice[time_index].grid.rho_tor_norm
        cs1d.grid.vol = vol = eq.time_slice[time_index].grid.vol

        auxHeatingProfile =  exp.(-4.0 * rho)
        pow_prof = cumul_integrate(vol, auxHeatingProfile)
        pow_prof = pow_prof ./ pow_prof[end]

        cs1d.electrons.power_inside = pow_prof .* Paux_e
        cs1d.total_ion_power_inside = pow_prof .* Paux_i
    end
end

function init(cpt, eqt, neped, ne_peaking, Te_ped, Te_peaking, w_ped, zeff, n_points=101)
    cpt.grid.rho_tor_norm = rho =  LinRange(0, 1, n_points)

    cpt.zeff = ones(inputs.rgrid) .* zeff
    cpt.rotation_frequency_tor_sonic = 5e3 .* abs(inputs.Paux .* 1.0 + 0.5) .* (1.0 .- rho)

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

    # Set densities
    ne_core = ne_peaking * neped
    ne = TAUENN_AD.Hmode_profiles(0.5 * neped, neped, ne_core, length(rho), 0.9, 0.9, w_ped)
    prof1d.electrons.density = ne
    zimp1 = 6.0
    niFraction = zeros(2)
    niFraction[2] = (zeff - 1.0) / (zimp1 * (zimp1 - 1.0))
    niFraction[1] = 1.0 - zimp1 * niFraction[2]
    for i in 1:length(prof1d.ion)
        prof1d.ion[i].density = ni = ne .* niFraction[i]
    end

    # Set temperatures
    # Use same approach of EPED guess for Tcore
    betaN = eqt...
    Bt = eqt...
    Ip = eqt...
    a = eqt...

    tcore = 10. * betaN * abs(Bt * Ip) / a / ne_core / (2.0 * 1.6e1 * 4.0 * pi * 1.0e-4)
    Te = Hmode_profiles(80., teped, tcore, length(rho), Te_peaking, Te_peaking, w_ped)
    prof1d.electrons.temperature = Te

    # to be done as an expression
    # prof1d.pressure_thermal = 1.6e-19 .* ne .* Te
    # for i in 1:length(prof1d.ion)
    #     ni = prof1d.ion[i].density
    #     prof1d.ion[i].temperature = Ti = tval .* inputs.Tratio
    #     pion =  1.6e-19 .* ni .* Ti
    #     prof1d.pressure_thermal = prof1d.pressure_thermal .+ pion
    # end

end