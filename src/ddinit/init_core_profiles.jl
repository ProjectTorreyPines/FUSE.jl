function init_core_profiles(dd::IMAS.dd; par::Parameters)
    init_core_profiles(
        dd.core_profiles,
        dd.equilibrium,
        dd.summary;
        ne_ped = par.core_profiles.ne_ped,
        ne_peaking = par.core_profiles.ne_peaking,
        Te_ped = par.core_profiles.Te_ped,
        Te_peaking = par.core_profiles.Te_peaking,
        w_ped = par.core_profiles.w_ped,
        zeff = par.core_profiles.zeff,
        P_co_nbi = par.core_profiles.P_co_nbi,
        ngrid = par.core_profiles.ngrid,
        bulk = par.core_profiles.bulk,
        impurity = par.core_profiles.impurity)
    return dd
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
    P_co_nbi::Real,
    T_ratio::Real = 1.0,
    ngrid::Int = 101
)
    cpt = resize!(cp.profiles_1d)

    cpt.grid.rho_tor_norm = LinRange(0, 1, ngrid)
    cpt.zeff = ones(ngrid) .* zeff
    cpt.rotation_frequency_tor_sonic = 5e3 * abs(P_co_nbi / 1e6 * 1.0 + 0.5) .* (1.0 .- cpt.grid.rho_tor_norm)

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
    cpt.electrons.density = AD_TAUENN.Hmode_profiles(0.5 * ne_ped, ne_ped, ne_core, ngrid, 0.9, 0.9, w_ped)

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
    cpt.electrons.temperature = AD_TAUENN.Hmode_profiles(80.0, Te_ped, Te_core, ngrid, Te_peaking, Te_peaking, w_ped)
    for i = 1:length(cpt.ion)
        cpt.ion[i].temperature = cpt.electrons.temperature ./ T_ratio
    end

    return cp
end
