"""
    init_core_profiles!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.core_profiles` starting from `ini` and `act` parameters
"""
function init_core_profiles!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_core_profiles")
    TimerOutputs.@timeit timer "init_core_profiles" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if IMAS.hasdata(dd1.core_profiles, :time) && length(dd1.core_profiles.time) > 0
                dd.core_profiles = deepcopy(dd1.core_profiles)
                # also set the pedestal in summary IDS
                if any([ismissing(getproperty(dd1.summary.local.pedestal, field), :value) for field in (:n_e, :zeff, :t_e)])
                    pe_ped, w_ped = IMAS.pedestal_finder(dd.core_profiles.profiles_1d[].electrons.pressure, dd.core_profiles.profiles_1d[].grid.psi_norm)
                    ped_summ = dd.summary.local.pedestal
                    cp1d = dd.core_profiles.profiles_1d[]
                    @ddtime ped_summ.position.rho_tor_norm = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.grid.rho_tor_norm).(1 - w_ped)
                    if ismissing(getproperty(dd1.summary.local.pedestal.n_e, :value, missing))
                        @ddtime ped_summ.n_e.value =
                            IMAS.interp1d(dd.core_profiles.profiles_1d[].grid.rho_tor_norm, dd.core_profiles.profiles_1d[].electrons.density_thermal).(1 - w_ped)
                    end
                    if ismissing(getproperty(dd1.summary.local.pedestal.t_e, :value, missing))
                        @ddtime ped_summ.t_e.value =
                            IMAS.interp1d(dd.core_profiles.profiles_1d[].grid.rho_tor_norm, dd.core_profiles.profiles_1d[].electrons.temperature).(1 - w_ped)
                    end
                    if ismissing(getproperty(dd1.summary.local.pedestal.t_i_average, :value, missing))
                        @ddtime ped_summ.t_i_average.value = IMAS.interp1d(dd.core_profiles.profiles_1d[].grid.rho_tor_norm, dd.core_profiles.profiles_1d[].t_i_average).(1 - w_ped)
                    end
                    if ismissing(getproperty(dd1.summary.local.pedestal.zeff, :value, missing))
                        @ddtime ped_summ.zeff.value = IMAS.interp1d(dd.core_profiles.profiles_1d[].grid.rho_tor_norm, dd.core_profiles.profiles_1d[].zeff).(1 - w_ped)
                    end
                end
            else
                init_from = :scalars
            end
            if ismissing(dd.core_profiles.global_quantities, :ejima) && !ismissing(ini.core_profiles, :ejima)
                @ddtime(dd.core_profiles.global_quantities.ejima = ini.core_profiles.ejima)
            end
        end

        if init_from == :scalars

            # if ini.equilibrium.pressure_core is not set, then estimate from ini.requirements.power_electric_net
            if ismissing(ini.equilibrium, :pressure_core) && !ismissing(ini.requirements, :power_electric_net)
                Pfusion_estimate = ini.requirements.power_electric_net * 2.0
                res = Optim.optimize(x -> cost_Pfusion_p0(x, Pfusion_estimate, dd, ini), 1e1, 1e7, Optim.GoldenSection())
                ini.equilibrium.pressure_core = pressure_core = res.minimizer[1]
                ActorCurrent(dd, act; ip_from=:pulse_schedule)
                ActorEquilibrium(dd, act; ip_from=:core_profiles)
            else
                pressure_core = ini.equilibrium.pressure_core
            end

            init_core_profiles!(
                dd.core_profiles,
                dd.equilibrium,
                dd.summary;
                pressure_core,
                helium_fraction=ini.core_profiles.bulk == :D ? 0.0 : ini.core_profiles.helium_fraction,
                ini.core_profiles.ne_setting,
                ini.core_profiles.ne_value,
                ini.core_profiles.T_ratio,
                ini.core_profiles.T_shaping,
                ini.core_profiles.n_shaping,
                ini.core_profiles.ne_sep_to_ped_ratio,
                ini.core_profiles.w_ped,
                ini.core_profiles.zeff,
                ini.core_profiles.rot_core,
                ini.core_profiles.ngrid,
                ini.core_profiles.bulk,
                ini.core_profiles.impurity,
                ejima=getproperty(ini.core_profiles, :ejima, missing),
                ini.core_profiles.polarized_fuel_fraction)
        end

        return dd
    end
end

function init_core_profiles!(
    cp::IMAS.core_profiles,
    eq::IMAS.equilibrium,
    summary::IMAS.summary;
    ne_setting::Symbol,
    ne_value::Real,
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
    ne_sep_to_ped_ratio::Real = 0.25,
    ngrid::Int=101)

    cp1d = resize!(cp.profiles_1d)
    eqt = eq.time_slice[]

    cp1d.grid.rho_tor_norm = range(0, 1, ngrid)
    cp1d.zeff = ones(ngrid) .* zeff
    cp1d.rotation_frequency_tor_sonic = IMAS.Hmode_profiles(0.0, rot_core / 8, rot_core, length(cp1d.grid.rho_tor_norm), 1.4, 1.4, 0.05)

    # Density handling for H-mode profile
    if ne_setting == :ne_ped
        ne_ped = ne_value
        ne_core = 1.4 * ne_ped
    elseif ne_setting == :greenwald_fraction_ped
        ne_ped = ne_value * IMAS.greenwald_density(eqt)
        ne_core = 1.4 * ne_ped
    elseif ne_setting == :greenwald_fraction
        ne_ped = IMAS.greenwald_density(eqt) * 0.8 * ne_value
        function cost_greenwald_fraction(ne0, greenwald_fraction_wanted)
            cp1d.electrons.density_thermal = IMAS.Hmode_profiles(ne_sep_to_ped_ratio * ne_ped, ne_ped, ne0, ngrid, n_shaping, n_shaping, w_ped)
            return (IMAS.greenwald_fraction(eqt, cp1d) - greenwald_fraction_wanted)^2
        end
        res = Optim.optimize(x-> cost_greenwald_fraction(x,ne_value), ne_ped/10 ,ne_ped*10, Optim.GoldenSection(),rel_tol=1E-3)
        ne_core = res.minimizer 
    end
    # Set ions:
    bulk_ion, imp_ion, he_ion = resize!(cp1d.ion, 3)
    # 1. DT
    IMAS.ion_element!(bulk_ion, bulk)
    @assert bulk_ion.element[1].z_n == 1.0 "Bulk ion `$bulk` must be a Hydrogenic isotope [:H, :D, :DT, :T]"
    # 2. Impurity
    IMAS.ion_element!(imp_ion, impurity)
    # 3. He
    IMAS.ion_element!(he_ion, :He4)

    # pedestal
    @ddtime summary.local.pedestal.n_e.value = ne_ped
    @ddtime summary.local.pedestal.position.rho_tor_norm = 1 - w_ped
    @ddtime summary.local.pedestal.zeff.value = zeff

    # Set densities
    if ne_core < ne_ped
        @warn "The core density is lower than the pedestal density, lower the pedestal density (ini.core_profiles.ne_ped)"
    end
    cp1d.electrons.density_thermal = IMAS.Hmode_profiles(ne_sep_to_ped_ratio * ne_ped, ne_ped, ne_core, ngrid, n_shaping, n_shaping, w_ped)
    # Zeff and quasi neutrality for a helium constant fraction with one impurity specie
    niFraction = zeros(3)
    # DT == 1
    # Imp == 2
    # He == 3
    zimp = imp_ion.element[1].z_n
    niFraction[3] = helium_fraction
    niFraction[1] = (zimp - zeff + 4 * niFraction[3] - 2 * zimp * niFraction[3]) / (zimp - 1)
    niFraction[2] = (zeff - niFraction[1] - 4 * niFraction[3]) / zimp^2
    @assert !any(niFraction .< 0.0) "zeff impossible to match for given helium fraction [$helium_fraction] and zeff [$zeff]"
    ni_core = 0.0
    for i in 1:length(cp1d.ion)
        cp1d.ion[i].density_thermal = cp1d.electrons.density_thermal .* niFraction[i]
        ni_core += cp1d.electrons.density_thermal[1] * niFraction[i]
    end

    # Set temperatures
    Te_core = minimum([pressure_core / (ni_core + ne_core) / IMAS.constants.e, 20e3])
    Te_ped = sqrt(Te_core / 1000.0 / 3.0) * 1000.0
    @ddtime summary.local.pedestal.t_e.value = Te_ped

    cp1d.electrons.temperature = IMAS.Hmode_profiles(80.0, Te_ped, Te_core, ngrid, T_shaping, T_shaping, w_ped)
    for i in 1:length(cp1d.ion)
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


function cost_Pfusion_p0(p::Real, target_pfus::Real, dd::IMAS.dd, ini::ParametersAllInits)
    init_core_profiles!(
        dd.core_profiles,
        dd.equilibrium,
        dd.summary;
        greenwald_fraction=getproperty(ini.core_profiles, :greenwald_fraction, missing),
        greenwald_fraction_ped=getproperty(ini.core_profiles, :greenwald_fraction_ped, missing),
        ne_ped=getproperty(ini.core_profiles, :ne_ped, missing),
        pressure_core=p,
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
    return abs(IMAS.fusion_power(dd.core_profiles.profiles_1d[]) - target_pfus)
end