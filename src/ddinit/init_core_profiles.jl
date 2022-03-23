function init_core_profiles(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :gasc
        gasc = GASC(par.gasc.filename, par.gasc.case)
        init_core_profiles(dd, gasc; bulk = par.core_profiles.bulk)

    elseif init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        if !ismissing(dd1.core_profiles, :time) && length(keys(dd1.core_profiles.time)) > 0
            dd.global_time = max(dd.global_time, maximum(dd1.core_profiles.time))
            dd.core_profiles = dd1.core_profiles
        else
            init_from = :scalars
        end
    end

    if init_from == :scalars
        init_core_profiles(
            dd.core_profiles,
            dd.equilibrium,
            dd.summary;
            ne_ped = par.core_profiles.ne_ped,
            n_peaking = par.core_profiles.n_peaking,
            T_shaping = par.core_profiles.T_shaping,
            w_ped = par.core_profiles.w_ped,
            zeff = par.core_profiles.zeff,
            rot_core = par.core_profiles.rot_core,
            ngrid = par.core_profiles.ngrid,
            bulk = par.core_profiles.bulk,
            impurity = par.core_profiles.impurity)
        @ddtime dd.core_profiles.global_quantities.ejima = par.core_profiles.ejima
    end

    return dd
end

function init_core_profiles(dd::IMAS.dd, gasc::GASC; bulk = :DT)
    gascsol = gasc.solution

    cp = dd.core_profiles
    cpt = resize!(cp.profiles_1d)

    cpt.grid.rho_tor_norm = gascsol["OUTPUTS"]["numerical profiles"]["rProf"]
    cpt.zeff = gascsol["OUTPUTS"]["numerical profiles"]["ZeffProf"]
    cpt.rotation_frequency_tor_sonic = cpt.grid.rho_tor_norm * 0.0 # < GASC has no notion of rotation

    # Set ions
    ion = resize!(cpt.ion, "label" => String(bulk))
    fill!(ion, IMAS.ion_element(bulk))
    @assert ion.element[1].z_n == 1 "Bulk ion must be a Hydrogen isotope [:H, :D, :DT, :T]"
    ion = resize!(cpt.ion, 2)
    element = resize!(ion.element, 1)
    element.z_n = gascsol["INPUTS"]["impurities"]["impurityZ"]
    element.a = Int(ceil(gascsol["INPUTS"]["impurities"]["impurityZ"] * 2.0))
    ion.label = "Impurity"

    # pedestal
    @ddtime dd.summary.local.pedestal.n_e.value = gascsol["OUTPUTS"]["plasma parameters"]["neped"] * 1E20
    i_ped = argmin(abs.(gascsol["OUTPUTS"]["numerical profiles"]["neProf"] .- gascsol["OUTPUTS"]["plasma parameters"]["neped"] / gascsol["OUTPUTS"]["plasma parameters"]["ne0"]))
    rho_ped = gascsol["OUTPUTS"]["numerical profiles"]["rProf"][i_ped]
    @ddtime dd.summary.local.pedestal.position.rho_tor_norm = rho_ped
    @ddtime dd.summary.local.pedestal.zeff.value = cpt.zeff[i_ped]

    # Set densities
    cpt.electrons.density = gascsol["OUTPUTS"]["numerical profiles"]["neProf"] * gascsol["OUTPUTS"]["plasma parameters"]["ne0"] * 1E20
    zimp1 = gascsol["INPUTS"]["impurities"]["impurityZ"]
    niFraction = zeros(2, length(cpt.grid.rho_tor_norm))
    niFraction[2, :] .= (cpt.zeff .- 1.0) ./ (zimp1 * (zimp1 - 1.0))
    niFraction[1, :] .= 1.0 .- zimp1 .* niFraction[2]
    @assert all(niFraction .> 0.0) "zeff too high for the impurity species with Z=" * string(gascsol["INPUTS"]["impurities"]["impurityZ"])
    for i = 1:length(cpt.ion)
        cpt.ion[i].density = cpt.electrons.density .* niFraction[i]
    end

    # Set temperatures
    Ti = gascsol["OUTPUTS"]["numerical profiles"]["TiProf"] * gascsol["INPUTS"]["plasma parameters"]["Ti0"] * 1E3
    cpt.electrons.temperature = Ti * gascsol["INPUTS"]["plasma parameters"]["Tratio"]
    for i = 1:length(cpt.ion)
        cpt.ion[i].temperature = Ti
    end

    # ejima
    IMAS.set_time_array(dd.core_profiles.global_quantities, :ejima, gascsol["INPUTS"]["plasma parameters"]["ejimaCoeff"])

    return dd
end

function init_core_profiles(
    cp::IMAS.core_profiles,
    eq::IMAS.equilibrium,
    summary::IMAS.summary;
    ne_ped::Real,
    n_peaking::Real,
    w_ped::Real,
    zeff::Real,
    bulk::Symbol,
    impurity::Symbol,
    rot_core::Real,
    T_ratio::Real = 1.0,
    T_shaping::Real = 1.8,
    n_shaping::Real = 0.9,
    ngrid::Int = 101
)
    cpt = resize!(cp.profiles_1d)

    cpt.grid.rho_tor_norm = LinRange(0, 1, ngrid)
    cpt.zeff = ones(ngrid) .* zeff
    cpt.rotation_frequency_tor_sonic = rot_core .* (1.0 .- cpt.grid.rho_tor_norm)

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
    ne_core = n_peaking * ne_ped
    cpt.electrons.density = TAUENN.Hmode_profiles(0.5 * ne_ped, ne_ped, ne_core, ngrid, n_shaping, n_shaping, w_ped)

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
    Te_ped = Te_core / 4
    cpt.electrons.temperature = TAUENN.Hmode_profiles(80.0, Te_ped, Te_core, ngrid, T_shaping, T_shaping, w_ped)
    for i = 1:length(cpt.ion)
        cpt.ion[i].temperature = cpt.electrons.temperature ./ T_ratio
    end

    return cp
end
