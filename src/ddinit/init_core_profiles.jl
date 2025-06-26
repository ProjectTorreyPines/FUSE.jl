"""
    init_core_profiles!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.core_profiles` starting from `ini` and `act` parameters
"""
function init_core_profiles!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_core_profiles")
    TimerOutputs.@timeit timer "init_core_profiles" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if !isempty(dd1.core_profiles.profiles_1d)
                dd.core_profiles = deepcopy(dd1.core_profiles)
                cp1d = dd.core_profiles.profiles_1d[]

                # also set the pedestal in summary IDS
                if any([ismissing(getproperty(dd.summary.local.pedestal, field), :value) for field in (:n_e, :zeff, :t_e)])
                    pedestal = IMAS.pedestal_finder(cp1d.electrons.pressure, cp1d.grid.psi_norm)
                    ped_summ = dd.summary.local.pedestal
                    @ddtime ped_summ.position.rho_tor_norm = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.grid.rho_tor_norm).(1 - pedestal.width)
                    if ismissing(getproperty(dd1.summary.local.pedestal.n_e, :value, missing))
                        @ddtime ped_summ.n_e.value =
                            IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(1 - pedestal.width)
                    end
                    if ismissing(getproperty(dd1.summary.local.pedestal.t_e, :value, missing))
                        @ddtime ped_summ.t_e.value =
                            IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(1 - pedestal.width)
                    end
                    if ismissing(getproperty(dd1.summary.local.pedestal.t_i_average, :value, missing))
                        @ddtime ped_summ.t_i_average.value = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average).(1 - pedestal.width)
                    end
                    if ismissing(getproperty(dd1.summary.local.pedestal.zeff, :value, missing))
                        @ddtime ped_summ.zeff.value = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.zeff).(1 - pedestal.width)
                    end
                end
            else
                init_from = :scalars
            end

            # ejima
            if ismissing(dd.core_profiles.global_quantities, :ejima) && !ismissing(ini.core_profiles, :ejima)
                @ddtime(dd.core_profiles.global_quantities.ejima = ini.core_profiles.ejima)
            end
        end

        if ini.core_profiles.ne_setting in (:ne_ped, :greenwald_fraction_ped)
            @assert act.ActorPedestal.density_match == :ne_ped "ini.core_profiles.ne_setting=:$(ini.core_profiles.ne_setting) requires act.ActorPedestal.density_match=:ne_ped"
        else
            @assert act.ActorPedestal.density_match == :ne_line "ini.core_profiles.ne_setting=:$(ini.core_profiles.ne_setting) requires act.ActorPedestal.density_match=:ne_line"
        end

        if !isempty(dd.equilibrium.time_slice)
            equil = dd.equilibrium.time_slice[]
        else
            equil = ini.equilibrium
        end

        if init_from == :scalars
            init_core_profiles!(
                dd.core_profiles,
                equil,
                dd.summary;
                ini.core_profiles.plasma_mode,
                ini.core_profiles.ne_setting,
                ini.core_profiles.ne_value,
                ini.core_profiles.ne_sep_to_ped_ratio,
                ini.core_profiles.ne_core_to_ped_ratio,
                ini.core_profiles.ne_shaping,
                pressure_core=getproperty(ini.equilibrium, :pressure_core, missing),
                helium_fraction=ini.core_profiles.bulk == :D ? 0.0 : ini.core_profiles.helium_fraction,
                ini.core_profiles.w_ped,
                ini.core_profiles.zeff,
                ini.core_profiles.bulk,
                ini.core_profiles.impurity,
                ini.core_profiles.rot_core,
                ejima=getproperty(ini.core_profiles, :ejima, missing),
                ini.core_profiles.polarized_fuel_fraction,
                ini.core_profiles.Ti_Te_ratio,
                ini.core_profiles.Te_shaping,
                ini.core_profiles.Te_sep,
                Te_ped=getproperty(ini.core_profiles, :Te_ped, missing),
                Te_core=getproperty(ini.core_profiles, :Te_core, missing),
                ini.core_profiles.ngrid,
                ITB_radius=getproperty(ini.core_profiles.ITB, :radius, missing),
                ITB_ne_width=getproperty(ini.core_profiles.ITB, :ne_width, missing),
                ITB_ne_height_ratio=getproperty(ini.core_profiles.ITB, :ne_height_ratio, missing),
                ITB_Te_width=getproperty(ini.core_profiles.ITB, :Te_width, missing),
                ITB_Te_height_ratio=getproperty(ini.core_profiles.ITB, :Te_height_ratio, missing)
            )
        end

        return dd
    end
end

function cost_WPED_α(rho::AbstractVector{<:Real}, profile0::AbstractVector{<:Real}, α::Real, value::Real, rho_ped::Real)
    profile = deepcopy(profile0)
    return cost_WPED_α!(rho, profile, α, value, rho_ped)
end

function init_core_profiles!(
    cp::IMAS.core_profiles,
    equil::Union{IMAS.equilibrium__time_slice,FUSEparameters__equilibrium},
    summary::IMAS.summary;
    plasma_mode::Symbol,
    ne_setting::Symbol,
    ne_value::Real,
    ne_sep_to_ped_ratio::Real,
    ne_core_to_ped_ratio::Real,
    ne_shaping::Real,
    pressure_core::Union{Real,Missing},
    helium_fraction::Real,
    w_ped::Real,
    zeff::Real,
    bulk::Symbol,
    impurity::Symbol,
    rot_core::Real,
    ejima::Union{Real,Missing},
    polarized_fuel_fraction::Real,
    Ti_Te_ratio::Real,
    Te_shaping::Real,
    Te_sep::Real,
    Te_ped::Union{Real,Missing},
    Te_core::Union{Real,Missing},
    ngrid::Int,
    ITB_radius::Union{Real,Missing},
    ITB_ne_width::Union{Real,Missing},
    ITB_ne_height_ratio::Union{Real,Missing},
    ITB_Te_width::Union{Real,Missing},
    ITB_Te_height_ratio::Union{Real,Missing})

    @assert 0.0 < ne_sep_to_ped_ratio < 1.0
    @assert ne_core_to_ped_ratio > 1.0
    @assert plasma_mode in (:H_mode, :L_mode)

    cp1d = resize!(cp.profiles_1d)

    # grid
    cp1d.grid.rho_tor_norm = range(0, 1, ngrid)

    # Density
    # first we start with a "unit" density profile...
    # NOTE: Throughout FUSE, the "pedestal" density is the density at rho=0.9
    cp1d.electrons.density_thermal = IMAS.Hmode_profiles(ne_sep_to_ped_ratio, 1.0, ne_core_to_ped_ratio, ngrid, ne_shaping, ne_shaping, w_ped)
    cp1d.electrons.density_thermal = IMAS.ped_height_at_09(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, 1.0)
    # ...which then we scale according to :ne_setting and :ne_value
    if ne_setting == :ne_ped
        ne_ped = ne_value
    elseif ne_setting == :greenwald_fraction_ped
        if typeof(equil) <: IMAS.equilibrium__time_slice
            ne_ped = ne_value * IMAS.greenwald_density(equil)
        else
            ne_ped = ne_value * IMAS.greenwald_density(equil.ip, equil.ϵ * equil.R0)
        end
    elseif ne_setting == :ne_line
        if typeof(equil) <: IMAS.equilibrium__time_slice
            ne_ped = ne_value / IMAS.ne_line(equil, cp1d)
        else
            ne_ped = ne_value / IMAS.ne_line(nothing, cp1d)
        end
    elseif ne_setting == :greenwald_fraction
        if typeof(equil) <: IMAS.equilibrium__time_slice
            ne_ped = ne_value * IMAS.greenwald_density(equil) / IMAS.ne_line(equil, cp1d)
        else
            ne_ped = ne_value * IMAS.greenwald_density(equil.ip, equil.ϵ * equil.R0) / IMAS.ne_line(nothing, cp1d)
        end
    end
    cp1d.electrons.density_thermal .*= ne_ped
    if ITB_radius !== missing
        cp1d.electrons.density_thermal = IMAS.ITB_profile(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, ITB_radius, ITB_ne_width, ITB_ne_height_ratio)
    end
    ne_core = cp1d.electrons.density_thermal[1]

    # Set ions:
    cp1d.zeff = ones(ngrid) .* zeff
    bulk_ion, imp_ion, he_ion = resize!(cp1d.ion, 3)
    # 1. DT
    IMAS.ion_element!(bulk_ion, bulk)
    @assert IMAS.is_hydrogenic(bulk_ion) "Bulk ion `$bulk` must be a Hydrogenic isotope [:H, :D, :DT, :T]"
    # 2. Impurity
    IMAS.ion_element!(imp_ion, impurity)
    # 3. He
    IMAS.ion_element!(he_ion, :He4)

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
    for i in eachindex(cp1d.ion)
        cp1d.ion[i].density_thermal = cp1d.electrons.density_thermal .* niFraction[i]
        ni_core += cp1d.electrons.density_thermal[1] * niFraction[i]
    end
    # remove He if not present
    if sum(niFraction[3]) == 0.0
        deleteat!(cp1d.ion, 3)
    end

    # temperatures
    @assert Te_core !== missing || pressure_core !== missing "Must specify either `ini.core_profiles.Te_core` or `ini.equilibrium.pressure_core`"
    if Te_core === missing
        Te_core = pressure_core / (Ti_Te_ratio * ni_core + ne_core) / IMAS.mks.e
    end
    if Te_ped === missing
        Te_ped = (Te_core - Te_sep) * w_ped .+ Te_sep
    end
    if plasma_mode == :H_mode
        cp1d.electrons.temperature = IMAS.Hmode_profiles(Te_sep, Te_ped, Te_core, ngrid, Te_shaping, Te_shaping, w_ped)
    else
        cp1d.electrons.temperature = IMAS.Lmode_profiles(Te_sep, Te_ped, Te_core, ngrid, Te_shaping, 1.0, w_ped)
    end
    if ITB_radius !== missing
        cp1d.electrons.temperature = IMAS.ITB_profile(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature, ITB_radius, ITB_Te_width, ITB_Te_height_ratio)
    end
    for i in eachindex(cp1d.ion)
        cp1d.ion[i].temperature = cp1d.electrons.temperature .* Ti_Te_ratio
    end

    # pedestal
    @ddtime summary.local.pedestal.n_e.value = ne_ped
    @ddtime summary.local.pedestal.position.rho_tor_norm = 1 - w_ped
    @ddtime summary.local.pedestal.zeff.value = zeff
    @ddtime summary.local.pedestal.t_e.value = Te_ped
    @ddtime summary.local.pedestal.t_i_average.value = Te_ped * Ti_Te_ratio

    # rotation
    if plasma_mode == :H_mode
        cp1d.rotation_frequency_tor_sonic = IMAS.Hmode_profiles(0.0, rot_core / ne_core_to_ped_ratio, rot_core, length(cp1d.grid.rho_tor_norm), Te_shaping, 1.0, w_ped)
    else
        cp1d.rotation_frequency_tor_sonic = IMAS.Lmode_profiles(0.0, rot_core / ne_core_to_ped_ratio, rot_core, length(cp1d.grid.rho_tor_norm), Te_shaping, 1.0, w_ped)
    end

    # ejima
    if ejima !== missing
        IMAS.set_time_array(cp.global_quantities, :ejima, ejima)
    end

    # set spin polarization
    cp.global_quantities.polarized_fuel_fraction = polarized_fuel_fraction
    return cp
end
