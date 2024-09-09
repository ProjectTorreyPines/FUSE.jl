"""
    ini_from_ods!(ini::ParametersAllInits; restore_expressions::Bool)::IMAS.dd

The purpose of this function is to setting `ini` values based on what is in the ods
thus simplifying the logic of the init functions after it which only have to look at ini values
"""
function ini_from_ods!(ini::ParametersAllInits; restore_expressions::Bool)::IMAS.dd
    if ini.general.init_from != :ods
        # don't do anything if to ini and return an empty dd
        dd1 = IMAS.dd()

    else
        # ini.general.dd takes priority
        if !ismissing(ini.general, :dd)
            dd1 = ini.general.dd
        else
            dd1 = load_ods(ini)
        end

        # equilibrium
        if !isempty(dd1.equilibrium.time_slice)
            eqt = dd1.equilibrium.time_slice[]
            fw = IMAS.first_wall(dd1.wall)
            IMAS.flux_surfaces(eqt, fw.r, fw.z)
            if ismissing(ini.equilibrium, :R0) && !ismissing(dd1.equilibrium.vacuum_toroidal_field, :r0)
                ini.equilibrium.R0 = dd1.equilibrium.vacuum_toroidal_field.r0
            end
            if ismissing(ini.equilibrium, :B0) && !ismissing(dd1.equilibrium.vacuum_toroidal_field, :b0)
                ini.equilibrium.B0 = @ddtime dd1.equilibrium.vacuum_toroidal_field.b0
            end
            if ismissing(ini.equilibrium, :ϵ) && !ismissing(dd1.equilibrium.time_slice[].boundary, :minor_radius)
                ini.equilibrium.ϵ = dd1.equilibrium.time_slice[].boundary.minor_radius / ini.equilibrium.R0
            end
            if ismissing(ini.equilibrium, :κ) && !ismissing(dd1.equilibrium.time_slice[].boundary, :elongation)
                ini.equilibrium.κ = dd1.equilibrium.time_slice[].boundary.elongation
            end
            if ismissing(ini.equilibrium, :δ) && !ismissing(dd1.equilibrium.time_slice[].boundary, :triangularity)
                ini.equilibrium.δ = dd1.equilibrium.time_slice[].boundary.triangularity
            end
            if ismissing(ini.equilibrium, :pressure_core) && !ismissing(eqt.profiles_1d, :pressure)
                ini.equilibrium.pressure_core = eqt.profiles_1d.pressure[1]
            end
            if ismissing(ini.equilibrium, :ip) && !ismissing(eqt.global_quantities, :ip)
                ini.equilibrium.ip = eqt.global_quantities.ip
            end
            if ismissing(ini.equilibrium, :xpoints)
                # look for x-points that fall within the first wall (if first-wall info is available)
                x_points = IMAS.x_points_inside_wall(eqt.boundary.x_point, dd1.wall)
                upper = any(x_point.z > eqt.boundary.geometric_axis.z for x_point in x_points)
                lower = any(x_point.z < eqt.boundary.geometric_axis.z for x_point in x_points)
                if upper && lower
                    ini.equilibrium.xpoints = :double
                elseif upper && !lower
                    ini.equilibrium.xpoints = :upper
                elseif !upper && lower
                    ini.equilibrium.xpoints = :lower
                else
                    ini.equilibrium.xpoints = :none
                end
            end
        end

        # core_profiles
        if ismissing(ini.core_profiles, :ejima) && !ismissing(dd1.core_profiles.global_quantities, :ejima)
            ini.core_profiles.ejima = @ddtime(dd1.core_profiles.global_quantities.ejima)
        end
        if !isempty(dd1.core_profiles.profiles_1d) && !ismissing(dd1.core_profiles.profiles_1d[].electrons, :pressure)
            cp1d = dd1.core_profiles.profiles_1d[]
            pe_ped, w_ped = IMAS.pedestal_finder(cp1d.electrons.pressure, cp1d.grid.psi_norm)
            ne_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(1 - w_ped)
            te_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(1 - w_ped)
            ti_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average).(1 - w_ped)
            zeff_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.zeff).(1 - w_ped)

            if ismissing(ini.core_profiles, :w_ped)
                ini.core_profiles.w_ped = w_ped
            end
            if ismissing(ini.core_profiles, :zeff)
                ini.core_profiles.zeff = zeff_ped
            end
            if ismissing(ini.core_profiles, :T_ratio)
                ini.core_profiles.T_ratio = ti_ped / te_ped
            end
        end

        # Here we delete fields from the ODS for which we know FUSE has expressions for.
        # Besides ensuring consistency, this is done because some FUSE workflows in fact expect certain fields to be expressions!
        if restore_expressions
            verbose = any(!contains(filename, "__FUSE__") for filename in split(ini.ods.filename, ","))
            FUSE.restore_init_expressions!(dd1; verbose)
        end

    end

    return dd1
end
