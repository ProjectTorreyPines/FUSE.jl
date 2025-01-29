"""
    ini_from_ods!(ini::ParametersAllInits; restore_expressions::Bool)::IMAS.dd

The purpose of this function is to set `ini` values based on what is in the ods thus
simplifying the logic of the init functions so that they only have to look at `ini` scalar values.
"""
function ini_from_ods!(ini::ParametersAllInits; restore_expressions::Bool)::IMAS.dd
    if ini.general.init_from != :ods
        # don't do anything if to ini and return an empty dd
        dd1 = IMAS.dd()

    else
        # note that ini.general.dd takes priority over ini.general.ods
        dd1 = load_ods(ini)

        # equilibrium
        eqt = nothing
        if !isempty(dd1.equilibrium.time_slice)
            eqt = dd1.equilibrium.time_slice[]

            if ismissing(ini.equilibrium, :B0) && !ismissing(dd1.equilibrium.vacuum_toroidal_field, :b0)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.B0 = @ddtime(dd1.equilibrium.vacuum_toroidal_field.b0)
                else
                    ini.equilibrium.B0 = t -> IMAS.moving_average(dd1.equilibrium.time, dd1.equilibrium.vacuum_toroidal_field.b0, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if ismissing(ini.equilibrium, :R0)
                if !ismissing(dd1.equilibrium.vacuum_toroidal_field, :r0)
                    ini.equilibrium.R0 = dd1.equilibrium.vacuum_toroidal_field.r0
                end
            end
            if ismissing(ini.equilibrium, :pressure_core) && !ismissing(eqt.profiles_1d, :pressure)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.pressure_core = eqt.profiles_1d.pressure[1]
                else
                    pressure_core = [eqt.profiles_1d.pressure[1] for eqt in dd1.equilibrium.time_slice]
                    ini.equilibrium.pressure_core = t -> IMAS.moving_average(dd1.equilibrium.time, pressure_core, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if ismissing(ini.equilibrium, :ip) && !ismissing(eqt.global_quantities, :ip)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.ip = eqt.global_quantities.ip
                else
                    ip = [eqt.global_quantities.ip for eqt in dd1.equilibrium.time_slice]
                    ini.equilibrium.ip = t -> IMAS.moving_average(dd1.equilibrium.time, ip, t, ini.time.pulse_shedule_time_basis.step)
                end
            end

            if !ismissing(eqt.boundary.outline, :r)
                mxh = IMAS.MXH(eqt.boundary.outline.r, eqt.boundary.outline.z)
                mxhs = [IMAS.MXH(eqt.boundary.outline.r, eqt.boundary.outline.z) for eqt in dd1.equilibrium.time_slice]
            else
                mxh = nothing
                mxhs = nothing
            end
            if ismissing(ini.equilibrium, :ϵ)
                if mxhs !== nothing
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.equilibrium.ϵ = mxh.ϵ
                    else
                        ϵ = [mxh.ϵ for mxh in mxhs]
                        ini.equilibrium.ϵ = t -> IMAS.moving_average(dd1.equilibrium.time, ϵ, t, ini.time.pulse_shedule_time_basis.step)
                    end
                elseif !ismissing(eqt.boundary, :minor_radius)
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.equilibrium.ϵ = eqt.boundary.minor_radius / ini.equilibrium.R0
                    else
                        data = [eqt.boundary.minor_radius / ini.equilibrium.R0 for eqt in dd1.equilibrium.time_slice]
                        ini.equilibrium.ϵ = t -> IMAS.moving_average(dd1.equilibrium.time, data, t, ini.time.pulse_shedule_time_basis.step)
                    end
                end
            end
            if ismissing(ini.equilibrium, :κ)
                if mxhs !== nothing
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.equilibrium.κ = mxh.κ
                    else
                        κ = [mxh.κ for mxh in mxhs]
                        ini.equilibrium.κ = t -> IMAS.moving_average(dd1.equilibrium.time, κ, t, ini.time.pulse_shedule_time_basis.step)
                    end
                elseif !ismissing(eqt.boundary, :elongation)
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.equilibrium.κ = eqt.boundary.elongation
                    else
                        elongation = [eqt.boundary.elongation for eqt in dd1.equilibrium.time_slice]
                        ini.equilibrium.κ = t -> IMAS.moving_average(dd1.equilibrium.time, elongation, t, ini.time.pulse_shedule_time_basis.step)
                    end
                end
            end
            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.tilt = mxh.tilt
                else
                    tilt = [mxh.tilt for mxh in mxhs]
                    ini.equilibrium.tilt = t -> IMAS.moving_average(dd1.equilibrium.time, tilt, t, ini.time.pulse_shedule_time_basis.step)
                end
            elseif !ismissing(eqt.boundary, :tilt)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.tilt = eqt.boundary.tilt
                else
                    tilt = [eqt.boundary.tilt for eqt in dd1.equilibrium.time_slice]
                    ini.equilibrium.tilt = t -> IMAS.moving_average(dd1.equilibrium.time, tilt, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.δ = mxh.δ
                else
                    δ = [mxh.δ for mxh in mxhs]
                    ini.equilibrium.δ = t -> IMAS.moving_average(dd1.equilibrium.time, δ, t, ini.time.pulse_shedule_time_basis.step)
                end
            elseif !ismissing(eqt.boundary, :triangularity)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.δ = eqt.boundary.triangularity
                else
                    triangularity = [eqt.boundary.triangularity for eqt in dd1.equilibrium.time_slice]
                    ini.equilibrium.δ = t -> IMAS.moving_average(dd1.equilibrium.time, triangularity, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.ζ = mxh.ζ
                else
                    ζ = [mxh.ζ for mxh in mxhs]
                    ini.equilibrium.ζ = t -> IMAS.moving_average(dd1.equilibrium.time, ζ, t, ini.time.pulse_shedule_time_basis.step)
                end
            elseif !ismissing(eqt.boundary, :squareness)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.ζ = eqt.boundary.squareness
                else
                    squareness = [eqt.boundary.squareness for eqt in dd1.equilibrium.time_slice]
                    ini.equilibrium.ζ = t -> IMAS.moving_average(dd1.equilibrium.time, squareness, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.𝚶 = mxh.𝚶
                else
                    𝚶 = [mxh.𝚶 for mxh in mxhs]
                    ini.equilibrium.𝚶 = t -> IMAS.moving_average(dd1.equilibrium.time, 𝚶, t, ini.time.pulse_shedule_time_basis.step)
                end
            elseif !ismissing(eqt.boundary, :ovality)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.𝚶 = eqt.boundary.ovality
                else
                    ovality = [eqt.boundary.ovality for eqt in dd1.equilibrium.time_slice]
                    ini.equilibrium.𝚶 = t -> IMAS.moving_average(dd1.equilibrium.time, ovality, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.twist = mxh.twist
                else
                    twist = [mxh.twist for mxh in mxhs]
                    ini.equilibrium.twist = t -> IMAS.moving_average(dd1.equilibrium.time, twist, t, ini.time.pulse_shedule_time_basis.step)
                end
            elseif !ismissing(eqt.boundary, :twist)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.twist = eqt.boundary.twist
                else
                    twist = [eqt.boundary.twist for eqt in dd1.equilibrium.time_slice]
                    ini.equilibrium.twist = t -> IMAS.moving_average(dd1.equilibrium.time, twist, t, ini.time.pulse_shedule_time_basis.step)
                end
            end

            if ismissing(ini.equilibrium, :xpoints)
                # look for x-points that fall within the first wall (if first-wall info is available)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    x_points = IMAS.x_points_inside_wall(eqt.boundary.x_point, dd1.wall)
                    upper = any(x_point.z > eqt.boundary.geometric_axis.z for x_point in x_points)
                    lower = any(x_point.z < eqt.boundary.geometric_axis.z for x_point in x_points)
                    ini.equilibrium.xpoints = xpoints_int_2_symbol(xpoints_bool_2_int(upper, lower))
                else
                    time_x_points = [IMAS.x_points_inside_wall(eqt.boundary.x_point, dd1.wall) for eqt in dd1.equilibrium.time_slice]
                    upper = [any(x_point.z > eqt.boundary.geometric_axis.z for x_point in x_points) for x_points in time_x_points]
                    lower = [any(x_point.z < eqt.boundary.geometric_axis.z for x_point in x_points) for x_points in time_x_points]
                    xpoints = Float64.(xpoints_bool_2_int.(upper, lower))
                    ini.equilibrium.xpoints = t -> xpoints_int_2_symbol(Int(round(IMAS.interp1d(dd1.equilibrium.time, xpoints).(t))))
                end
            end
        end

        # core_profiles
        if !isempty(dd1.core_profiles.profiles_1d)
            NE_PED = Float64[]
            W_PED = Float64[]
            ZEFF = Float64[]
            TI_TE_RATIO = Float64[]
            TE_CORE = Float64[]
            if ismissing(ini.time, :pulse_shedule_time_basis)
                time = [dd1.global_time]
            else
                time = dd1.core_profiles.time
                time = time[dd1.equilibrium.time[1].<=time.<=dd1.equilibrium.time[end]]
            end
            for time0 in time
                cp1d = dd1.core_profiles.profiles_1d[time0]
                eqt = dd1.equilibrium.time_slice[time0]
                if ismissing(cp1d.grid, :psi) && eqt !== nothing && !ismissing(eqt.profiles_1d, :rho_tor_norm)
                    cp1d.grid.psi = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.psi).(cp1d.grid.rho_tor_norm)
                end
                if !ismissing(cp1d.electrons, :pressure)
                    psi_norm = getproperty(cp1d.grid, :psi_norm, cp1d.grid.rho_tor_norm) # ideally `psi_norm`, but if not available `rho_tor_norm` will do
                    _, w_ped = IMAS.pedestal_finder(cp1d.electrons.pressure, psi_norm)
                    ne_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(1 - w_ped)
                    te_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(1 - w_ped)
                    ti_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average).(1 - w_ped)
                    zeff_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.zeff).(1 - w_ped)

                    push!(NE_PED, ne_ped)
                    push!(W_PED, w_ped)
                    push!(ZEFF, zeff_ped)
                    push!(TI_TE_RATIO, ti_ped / te_ped)
                end
                if !ismissing(cp1d.electrons, :temperature)
                    push!(TE_CORE, cp1d.electrons.temperature[1] / 1E3)
                end
            end
            if ismissing(ini.core_profiles, :ne_setting) && ismissing(ini.core_profiles, :ne_value)
                ini.core_profiles.ne_setting = :ne_ped
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.ne_value = NE_PED[1]
                else
                    ini.core_profiles.ne_value = t -> IMAS.moving_average(time, NE_PED, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if ismissing(ini.core_profiles, :w_ped)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.w_ped = W_PED[1]
                else
                    ini.core_profiles.w_ped = t -> IMAS.moving_average(time, W_PED, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if ismissing(ini.core_profiles, :zeff)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.zeff = ZEFF[1]
                else
                    ini.core_profiles.zeff = t -> IMAS.moving_average(time, ZEFF, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if ismissing(ini.core_profiles, :Ti_Te_ratio)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.Ti_Te_ratio = TI_TE_RATIO[1]
                else
                    ini.core_profiles.Ti_Te_ratio = t -> IMAS.moving_average(time, TI_TE_RATIO, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if ismissing(ini.core_profiles, :Te_core)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.Te_core = TE_CORE[1]
                else
                    ini.core_profiles.Te_core = t -> IMAS.moving_average(time, TE_CORE, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
            if ismissing(ini.core_profiles, :ejima) && !ismissing(dd1.core_profiles.global_quantities, :ejima)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.ejima = @ddtime dd1.core_profiles.global_quantities.ejima
                else
                    ini.core_profiles.ejima = t -> IMAS.moving_average(dd1.core_profiles.time, dd1.core_profiles.global_quantities.ejima, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
        end

        # EC
        for (k, beam) in enumerate(dd1.ec_launchers.beam)
            if k > length(ini.ec_launcher)
                resize!(ini.ec_launcher, k)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.ec_launcher[k].power_launched = @ddtime beam.power_launched.data
                else
                    ini.ec_launcher[k].power_launched = t -> IMAS.moving_average(beam.power_launched.time, beam.power_launched.data, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
        end

        # IC
        for (k, antenna) in enumerate(dd1.ic_antennas.antenna)
            if k > length(ini.ic_antenna)
                resize!(ini.ic_antenna, k)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.ic_antenna[k].power_launched = @ddtime antenna.power_launched.data
                else
                    ini.ic_antenna[k].power_launched = t -> IMAS.moving_average(antenna.power_launched.time, antenna.power_launched.data, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
        end

        # LH
        for (k, antenna) in enumerate(dd1.lh_antennas.antenna)
            if k > length(ini.lh_antenna)
                resize!(ini.lh_antenna, k)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.lh_antenna[k].power_launched = @ddtime antenna.power_launched.data
                else
                    ini.lh_antenna[k].power_launched = t -> IMAS.moving_average(antenna.power_launched.time, antenna.power_launched.data, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
        end

        # NB
        for (k, unit) in enumerate(dd1.nbi.unit)
            if k > length(ini.nb_unit)
                resize!(ini.nb_unit, k)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.nb_unit[k].power_launched = @ddtime unit.power_launched.data
                else
                    ini.nb_unit[k].power_launched = t -> IMAS.moving_average(unit.power_launched.time, unit.power_launched.data, t, ini.time.pulse_shedule_time_basis.step)
                end
                # make beam energy constant
                ini.nb_unit[k].beam_energy = maximum(unit.energy.data)
                ini.nb_unit[k].toroidal_angle = unit.beamlets_group[1].angle
            end
        end

        # PL
        for (k, launcher) in enumerate(dd1.pellets.launcher)
            if k > length(ini.pellet_launcher)
                resize!(ini.pellet_launcher, k)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.pellet_launcher[k].frequency = @ddtime launcher.frequency.data
                else
                    ini.pellet_launcher[k].frequency = t -> IMAS.moving_average(launcher.frequency.time, launcher.frequency.data, t, ini.time.pulse_shedule_time_basis.step)
                end
            end
        end

        # Here we delete fields from the ODS for which we know FUSE has expressions for.
        # Besides ensuring consistency, this is done because some FUSE workflows in fact expect certain fields to be expressions!
        if restore_expressions
            verbose = !ismissing(ini.ods, :filename) && any(!contains(filename, "__FUSE__") for filename in split(ini.ods.filename, ","))
            FUSE.restore_init_expressions!(dd1; verbose)
        end
    end

    return dd1
end

function xpoints_int_2_symbol(xpoint::Int)
    if xpoint == 2
        return :double
    elseif xpoint == 1
        return :upper
    elseif xpoint == -1
        return :lower
    elseif xpoint == 0
        return :none
    end
end

function xpoints_bool_2_int(upper::Bool, lower::Bool)
    if upper && lower
        return 2
    elseif upper && !lower
        return 1
    elseif !upper && lower
        return -1
    else
        return 0
    end
end
