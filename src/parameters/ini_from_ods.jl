"""
    set_ini_act_from_ods!(ini::ParametersAllInits, act::ParametersAllActors)

The purpose of this function is to set `ini` values based on what is in the ods thus
simplifying the logic of the init functions so that they only have to look at `ini` scalar values.
"""
function set_ini_act_from_ods!(ini::ParametersAllInits, act::ParametersAllActors)
    if ini.general.init_from != :ods
        # don't do anything to ini and return an empty dd
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
                    ini.equilibrium.B0 = TimeData(dd1.equilibrium.time, dd1.equilibrium.vacuum_toroidal_field.b0)
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
                    ini.equilibrium.pressure_core = TimeData(dd1.equilibrium.time, [eqt.profiles_1d.pressure[1] for eqt in dd1.equilibrium.time_slice])
                end
            end
            if ismissing(ini.equilibrium, :ip) && !ismissing(eqt.global_quantities, :ip)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.ip = eqt.global_quantities.ip
                else
                    ini.equilibrium.ip = TimeData(dd1.equilibrium.time, [eqt.global_quantities.ip for eqt in dd1.equilibrium.time_slice])
                end
            end

            if !ismissing(eqt.boundary.outline, :r)
                pr, pz = IMAS.resample_plasma_boundary(eqt.boundary.outline.r, eqt.boundary.outline.z; n_points=100)
                mxh = IMAS.MXH(pr, pz)
                mxhs = [IMAS.MXH(pr, pz) for eqt in dd1.equilibrium.time_slice]
            else
                mxh = nothing
                mxhs = nothing
            end

            if ismissing(ini.equilibrium, :ϵ)
                if mxhs !== nothing
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.equilibrium.ϵ = mxh.ϵ
                    else
                        ini.equilibrium.ϵ = TimeData(dd1.equilibrium.time, [mxh.ϵ for mxh in mxhs])
                    end
                elseif !ismissing(eqt.boundary, :minor_radius)
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.equilibrium.ϵ = eqt.boundary.minor_radius / ini.equilibrium.R0
                    else
                        ini.equilibrium.ϵ = TimeData(dd1.equilibrium.time, [eqt.boundary.minor_radius / ini.equilibrium.R0 for eqt in dd1.equilibrium.time_slice])
                    end
                end
            end
            if ismissing(ini.equilibrium, :κ)
                if mxhs !== nothing
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.equilibrium.κ = mxh.κ
                    else
                        ini.equilibrium.κ = TimeData(dd1.equilibrium.time, [mxh.κ for mxh in mxhs])
                    end
                elseif !ismissing(eqt.boundary, :elongation)
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.equilibrium.κ = eqt.boundary.elongation
                    else
                        ini.equilibrium.κ = TimeData(dd1.equilibrium.time, [eqt.boundary.elongation for eqt in dd1.equilibrium.time_slice])
                    end
                end
            end

            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.tilt = mxh.tilt
                else
                    ini.equilibrium.tilt = TimeData(dd1.equilibrium.time, [mxh.tilt for mxh in mxhs])
                end
            elseif !ismissing(eqt.boundary, :tilt)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.tilt = eqt.boundary.tilt
                else
                    ini.equilibrium.tilt = TimeData(dd1.equilibrium.time, [eqt.boundary.tilt for eqt in dd1.equilibrium.time_slice])
                end
            end

            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.δ = mxh.δ
                else
                    ini.equilibrium.δ = TimeData(dd1.equilibrium.time, [mxh.δ for mxh in mxhs])
                end
            elseif !ismissing(eqt.boundary, :triangularity)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.δ = eqt.boundary.triangularity
                else
                    ini.equilibrium.δ = TimeData(dd1.equilibrium.time, [eqt.boundary.triangularity for eqt in dd1.equilibrium.time_slice])
                end
            end

            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.ζ = mxh.ζ
                else
                    ζ = [mxh.ζ for mxh in mxhs]
                    ini.equilibrium.ζ = TimeData(dd1.equilibrium.time, ζ)
                end
            elseif !ismissing(eqt.boundary, :squareness)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.ζ = eqt.boundary.squareness
                else
                    ini.equilibrium.ζ = TimeData(dd1.equilibrium.time, [eqt.boundary.squareness for eqt in dd1.equilibrium.time_slice])
                end
            end

            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.𝚶 = mxh.𝚶
                else
                    ini.equilibrium.𝚶 = TimeData(dd1.equilibrium.time, [mxh.𝚶 for mxh in mxhs])
                end
            elseif !ismissing(eqt.boundary, :ovality)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.𝚶 = eqt.boundary.ovality
                else
                    ini.equilibrium.𝚶 = TimeData(dd1.equilibrium.time, [eqt.boundary.ovality for eqt in dd1.equilibrium.time_slice])
                end
            end

            if mxhs !== nothing
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.twist = mxh.twist
                else
                    ini.equilibrium.twist = TimeData(dd1.equilibrium.time, [mxh.twist for mxh in mxhs])
                end
            elseif !ismissing(eqt.boundary, :twist)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.equilibrium.twist = eqt.boundary.twist
                else
                    ini.equilibrium.twist = TimeData(dd1.equilibrium.time, [eqt.boundary.twist for eqt in dd1.equilibrium.time_slice])
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
                    xpoints = xpoints_int_2_symbol.(Int.(Float64.(xpoints_bool_2_int.(upper, lower))))
                    ini.equilibrium.xpoints = TimeData(dd1.equilibrium.time, xpoints)
                end
            end
        end

        # core_profiles
        if !isempty(dd1.core_profiles.profiles_1d)
            W_PED = Float64[]
            NE_PED = Float64[]
            NE_LINE = Float64[]
            ZEFF = Float64[]
            TI_TE_RATIO = Float64[]
            TE_CORE = Float64[]
            if ismissing(ini.time, :pulse_shedule_time_basis)
                time = [dd1.global_time]
            else
                time = dd1.core_profiles.time
                time = time[dd1.equilibrium.time[1].<=time.<=dd1.equilibrium.time[end]]
            end
            pedestal = nothing
            for time0 in time
                cp1d = dd1.core_profiles.profiles_1d[time0]
                if ismissing(cp1d.grid, :psi) && !isempty(dd1.equilibrium.time_slice) && !ismissing(eqt.profiles_1d, :rho_tor_norm)
                    eqt = dd1.equilibrium.time_slice[time0]
                    cp1d.grid.psi = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.psi).(cp1d.grid.rho_tor_norm)
                end
                if !ismissing(cp1d.electrons, :pressure)
                    rho09 = 0.9 # in FUSE pedestal quantities are defined at rho=0.9
                    pedestal = IMAS.pedestal_finder(cp1d.electrons.pressure, cp1d.grid.rho_tor_norm; guess=pedestal)
                    ne_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(rho09)
                    te_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(rho09)
                    ti_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average).(rho09)
                    zeff_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.zeff).(rho09)

                    if isempty(dd1.equilibrium.time_slice)
                        nel = IMAS.ne_line(nothing, cp1d)
                    else
                        eqt = dd1.equilibrium.time_slice[time0]
                        nel = IMAS.ne_line(eqt, cp1d)
                    end

                    push!(W_PED, pedestal.width)
                    push!(NE_PED, ne_ped)
                    push!(NE_LINE, nel)
                    push!(ZEFF, zeff_ped)
                    push!(TI_TE_RATIO, ti_ped / te_ped)
                end
                if !ismissing(cp1d.electrons, :temperature)
                    push!(TE_CORE, cp1d.electrons.temperature[1] / 1E3)
                end
            end
            if ismissing(ini.core_profiles, :ne_value)
                if ismissing(ini.core_profiles, :ne_setting) || ini.core_profiles.ne_setting == :ne_ped
                    ini.core_profiles.ne_setting = :ne_ped
                    @assert !any(isnan.(NE_PED))
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.core_profiles.ne_value = NE_PED[1]
                    else
                        ini.core_profiles.ne_value = TimeData(time, NE_PED)
                    end
                elseif ini.core_profiles.ne_setting == :ne_line
                    @assert !any(isnan.(NE_LINE))
                    if ismissing(ini.time, :pulse_shedule_time_basis)
                        ini.core_profiles.ne_value = NE_LINE[1]
                    else
                        ini.core_profiles.ne_value = TimeData(time, NE_LINE)
                    end
                end
            end
            if ismissing(ini.core_profiles, :w_ped)
                @assert !any(isnan.(W_PED))
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.w_ped = W_PED[1]
                else
                    ini.core_profiles.w_ped = TimeData(time, W_PED)
                end
            end
            if ismissing(ini.core_profiles, :zeff)
                @assert !any(isnan.(ZEFF))
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.zeff = ZEFF[1]
                else
                    ini.core_profiles.zeff = TimeData(time, ZEFF)
                end
            end
            if ismissing(ini.core_profiles, :Ti_Te_ratio)
                @assert !any(isnan.(TI_TE_RATIO))
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.Ti_Te_ratio = TI_TE_RATIO[1]
                else
                    ini.core_profiles.Ti_Te_ratio = TimeData(time, TI_TE_RATIO)
                end
            end
            if ismissing(ini.core_profiles, :Te_core)
                @assert !any(isnan.(TE_CORE))
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.Te_core = TE_CORE[1]
                else
                    ini.core_profiles.Te_core = TimeData(time, TE_CORE)
                end
            end
            if ismissing(ini.core_profiles, :ejima) && !ismissing(dd1.core_profiles.global_quantities, :ejima)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.core_profiles.ejima = @ddtime dd1.core_profiles.global_quantities.ejima
                else
                    ini.core_profiles.ejima = TimeData(dd1.core_profiles.time, dd1.core_profiles.global_quantities.ejima)
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
                    ini.ec_launcher[k].power_launched = TimeData(beam.power_launched.time, beam.power_launched.data)
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
                    ini.ic_antenna[k].power_launched = TimeData(antenna.power_launched.time, antenna.power_launched.data)
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
                    ini.lh_antenna[k].power_launched = TimeData(antenna.power_launched.time, antenna.power_launched.data)
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
                    ini.nb_unit[k].power_launched = TimeData(unit.power_launched.time, unit.power_launched.data)
                end
                # make beam energy constant
                ini.nb_unit[k].beam_energy = maximum(unit.energy.data)
                ini.nb_unit[k].normalized_tangency_radius = unit.beamlets_group[1].tangency_radius / (0.5 * (eqt.profiles_1d.r_inboard[end] + eqt.profiles_1d.r_outboard[end]))
                ini.nb_unit[k].offaxis = abs(unit.beamlets_group[1].angle) > (0.1 * deg)
            end
        end

        # PL
        for (k, launcher) in enumerate(dd1.pellets.launcher)
            if k > length(ini.pellet_launcher)
                resize!(ini.pellet_launcher, k)
                if ismissing(ini.time, :pulse_shedule_time_basis)
                    ini.pellet_launcher[k].frequency = @ddtime launcher.frequency.data
                else
                    ini.pellet_launcher[k].frequency = TimeData(launcher.frequency.time, launcher.frequency.data)
                end
            end
        end
    end

    consistent_ini_act!(ini, act)

    return dd1
end

"""
    consistent_ini_act!(ini::ParametersAllInits, act::ParametersAllActors)

Makes `ini` and `act` self-consistent and consistent with one another
"""
function consistent_ini_act!(ini::ParametersAllInits, act::ParametersAllActors)
    if !isempty(ini.ec_launcher)
        if isempty(act.ActorSimpleEC.actuator)
            resize!(act.ActorSimpleEC.actuator, length(ini.ec_launcher))
        else
            @assert length(act.ActorSimpleEC.actuator) == length(ini.ec_launcher) "length(act.ActorSimpleEC.actuator) = $(length(act.ActorSimpleEC.actuator)) must be equal to length(ini.ec_launcher)=$(length(ini.ec_launcher))"
        end
    end

    if !isempty(ini.ic_antenna)
        if isempty(act.ActorSimpleIC.actuator)
            resize!(act.ActorSimpleIC.actuator, length(ini.ic_antenna))
        else
            @assert length(act.ActorSimpleIC.actuator) == length(ini.ic_antenna) "length(act.ActorSimpleIC.actuator) = $(length(act.ActorSimpleIC.actuator)) must be equal to length(ini.ic_antenna)=$(length(ini.ic_antenna))"
        end
    end

    if !isempty(ini.lh_antenna)
        if isempty(act.ActorSimpleLH.actuator)
            resize!(act.ActorSimpleLH.actuator, length(ini.lh_antenna))
        else
            @assert length(act.ActorSimpleLH.actuator) == length(ini.lh_antenna) "length(act.ActorSimpleLH.actuator) = $(length(act.ActorSimpleLH.actuator)) must be equal to length(ini.lh_antenna)=$(length(ini.lh_antenna))"
        end
    end

    if !isempty(ini.nb_unit)
        if isempty(act.ActorSimpleNB.actuator)
            resize!(act.ActorSimpleNB.actuator, length(ini.nb_unit))
        else
            @assert length(act.ActorSimpleNB.actuator) == length(ini.nb_unit) "length(act.ActorSimpleNB.actuator) = $(length(act.ActorSimpleNB.actuator)) must be equal to length(ini.nb_unit)=$(length(ini.nb_unit))"
        end
    end

    if !isempty(ini.pellet_launcher)
        if isempty(act.ActorSimplePL.actuator)
            resize!(act.ActorSimplePL.actuator, length(ini.pellet_launcher))
        else
            @assert length(act.ActorSimplePL.actuator) == length(ini.pellet_launcher) "length(act.ActorSimplePL.actuator) = $(length(act.ActorSimplePL.actuator)) must be equal to length(ini.pellet_launcher)=$(length(ini.pellet_launcher))"
        end
    end
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
