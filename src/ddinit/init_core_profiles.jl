"""
    init_core_profiles(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.core_profiles` starting from `ini` and `act` parameters
"""
function init_core_profiles(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    TimerOutputs.reset_timer!("init_core_profiles")
    TimerOutputs.@timeit timer "init_core_profiles" begin
        init_from = ini.general.init_from

        if init_from == :ods
            dd1 = IMAS.json2imas(ini.ods.filename)
            if !ismissing(dd1.core_profiles, :time) && length(keys(dd1.core_profiles.time)) > 0
                dd.global_time = max(dd.global_time, maximum(dd1.core_profiles.time))
                dd.core_profiles = dd1.core_profiles
            else
                init_from = :scalars
            end
            if ismissing(dd.core_profiles.global_quantities, :ejima) && ~ismissing(ini.core_profiles, :ejima)
                IMAS.set_time_array(dd.core_profiles.global_quantities, :ejima, ini.core_profiles.ejima)
            end
        end

        if init_from == :scalars
            init_core_profiles(
                dd.core_profiles,
                dd.equilibrium,
                dd.summary;
                greenwald_fraction=getproperty(ini.core_profiles, :greenwald_fraction, missing),
                greenwald_fraction_ped=getproperty(ini.core_profiles, :greenwald_fraction_ped, missing),
                ne_ped=getproperty(ini.core_profiles, :ne_ped, missing),
                pressure_core=dd.equilibrium.time_slice[].profiles_1d.pressure[1],
                helium_fraction=ini.core_profiles.helium_fraction,
                T_ratio=ini.core_profiles.T_ratio,
                T_shaping=ini.core_profiles.T_shaping,
                n_shaping=ini.core_profiles.n_shaping,
                w_ped=ini.core_profiles.w_ped,
                zeff=ini.core_profiles.zeff,
                rot_core=ini.core_profiles.rot_core,
                ngrid=ini.core_profiles.ngrid,
                bulk=ini.core_profiles.bulk,
                impurity=ini.core_profiles.impurity,
                ejima=getproperty(ini.core_profiles, :ejima, missing),
                polarized_fuel_fraction=ini.core_profiles.polarized_fuel_fraction)
        end

        return dd
    end
end

function init_core_profiles(
    cp::IMAS.core_profiles,
    eq::IMAS.equilibrium,
    summary::IMAS.summary;
    greenwald_fraction::Union{Real,Missing},
    greenwald_fraction_ped::Union{Real,Missing},
    ne_ped::Union{Real,Missing},
    pressure_core::Real,
    helium_fraction::Real,
    w_ped::Real,
    zeff::Real,
    bulk::Symbol,
    impurity::Symbol,
    rot_core::Real,
    ejima::Union{Real,Missing},
    polarized_fuel_fraction::Real,
    T_ratio::Real=1.0,
    T_shaping::Real=1.8,
    n_shaping::Real=0.9,
    ngrid::Int=101)

    cp1d = resize!(cp.profiles_1d)
    eqt = eq.time_slice[]

    cp1d.grid.rho_tor_norm = LinRange(0, 1, ngrid)
    cp1d.zeff = ones(ngrid) .* zeff
    cp1d.rotation_frequency_tor_sonic = rot_core .* (1.0 .- cp1d.grid.rho_tor_norm)

    # Density handling
    @assert !ismissing(ne_ped) || !ismissing(greenwald_fraction) || !ismissing(greenwald_fraction_ped) "At least one of ini.core_profiles ne_ped / greenwald_fraction_ped / greenwald_fraction must be set"
    @assert !ismissing(ne_ped) ⊻ !ismissing(greenwald_fraction_ped) "One and only one of ini.core_profiles.ne_ped or ini.core_profiles.greenwald_fraction_ped must be set"
    if ismissing(ne_ped)
        ne_ped = greenwald_fraction_ped * IMAS.greenwald_density(eqt)
    end
    if ismissing(greenwald_fraction)
        # guess greewald fraction from ne_ped
        ne0_guess = ne_ped * 1.4
        cp1d.electrons.density_thermal = IMAS.Hmode_profiles(0.5 * ne_ped, ne_ped, ne0_guess, ngrid, n_shaping, n_shaping, w_ped)
        greenwald_fraction = IMAS.greenwald_fraction(eqt, cp1d)
    end

    # Set ions:
    # 1. DT
    ion = resize!(cp1d.ion, "label" => String(bulk))
    fill!(ion, IMAS.ion_element(ion_symbol=bulk))
    @assert ion.element[1].z_n == 1.0 "Bulk ion must be a Hydrogenic isotope [:H, :D, :DT, :T]"

    # 2. Impurity
    ion = resize!(cp1d.ion, "label" => String(impurity))
    fill!(ion, IMAS.ion_element(ion_symbol=impurity))

    # 3. He
    ion = resize!(cp1d.ion, "label" => "He")
    fill!(ion, IMAS.ion_element(ion_symbol=:He))

    # pedestal
    @ddtime summary.local.pedestal.n_e.value = ne_ped
    @ddtime summary.local.pedestal.position.rho_tor_norm = 1 - w_ped
    @ddtime summary.local.pedestal.zeff.value = zeff

    # Set densities
    function cost_greenwald_fraction(ne0)
        ne0 = ne0[1]
        cp1d.electrons.density_thermal = IMAS.Hmode_profiles(0.5 * ne_ped, ne_ped, ne0, ngrid, n_shaping, n_shaping, w_ped)
        return (IMAS.greenwald_fraction(eqt, cp1d) - greenwald_fraction)^2
    end
    ne0_guess = ne_ped * 1.4
    res = Optim.optimize(cost_greenwald_fraction, [ne0_guess], Optim.NelderMead(), Optim.Options(g_tol=1E-4))
    ne_core = res.minimizer[1]
    if ne_core < ne_ped
        @warn "The core density is lower than the pedestal density, lower the pedestal density (ini.core_profiles.ne_ped)"
    end
    cp1d.electrons.density_thermal = IMAS.Hmode_profiles(0.5 * ne_ped, ne_ped, ne_core, ngrid, n_shaping, n_shaping, w_ped)
    # Zeff and quasi neutrality for a helium constant fraction with one impurity specie
    niFraction = zeros(3)
    # DT == 1
    # Imp == 2
    # He == 3
    zimp = IMAS.ion_element(ion_symbol=impurity).element[1].z_n
    niFraction[3] = helium_fraction
    niFraction[1] = (zimp - zeff + 4 * niFraction[3] - 2 * zimp * niFraction[3]) / (zimp - 1)
    niFraction[2] = (zeff - niFraction[1] - 4 * niFraction[3]) / zimp^2
    @assert !any(niFraction .< 0.0) "zeff impossible to match for given helium fraction [$helium_fraction] and zeff [$zeff]"
    ni_core = 0.0
    for i = 1:length(cp1d.ion)
        cp1d.ion[i].density_thermal = cp1d.electrons.density_thermal .* niFraction[i]
        ni_core += cp1d.electrons.density_thermal[1] * niFraction[i]
    end

    # Set temperatures
    Te_core = pressure_core / (ni_core + ne_core) / IMAS.constants.e
    Te_ped = sqrt(Te_core / 1000.0 / 3.0) * 1000.0
    @ddtime summary.local.pedestal.t_e.value = Te_ped

    cp1d.electrons.temperature = IMAS.Hmode_profiles(80.0, Te_ped, Te_core, ngrid, T_shaping, T_shaping, w_ped)
    for i = 1:length(cp1d.ion)
        cp1d.ion[i].temperature = cp1d.electrons.temperature .* T_ratio
    end

    # remove He if not present
    if sum(niFraction[3]) == 0.0
        deleteat!(cp1d.ion, 3)
    end

    # ejima
    if ejima !== missing
        IMAS.set_time_array(cp.global_quantities, :ejima, ejima)
    end

    # set spin polarization
    cp.global_quantities.polarized_fuel_fraction = polarized_fuel_fraction
    return cp
end
