"""
    ini_from_ods!(ini::ParametersAllInits)::IMAS.dd

The purpose of this function is to setting `ini` values based on what is in the ods
thus simplifying the logic of the init functions after it which only have to look at ini values
"""
function ini_from_ods!(ini::ParametersAllInits)::IMAS.dd
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
            IMAS.flux_surfaces(eqt)
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
                x_points = IMAS.x_points_in_wall(eqt.boundary.x_point, dd1.wall)
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
        if !ismissing(dd1.core_profiles.global_quantities, :ejima)
            ini.core_profiles.ejima = @ddtime(dd1.core_profiles.global_quantities.ejima)
        end

        # here we selectively copy the quantities that FUSE considers `primary`
        # all other quantities are either of no interest for intialization or
        # they have expressions associated with them.
        dd1 = copy_init_primary_quanties(dd1)

    end

    return dd1
end

"""
    copy_init_primary_quanties(dd::IMAS.DD)

here we selectively copy the quantities that FUSE considers `primary` to a new dd
All other quantities are either of no interest for intialization or they have expressions associated with them
"""
function copy_init_primary_quanties(dd::IMAS.DD)
    dd_out = typeof(dd)()
    dd_out.global_time = dd.global_time

    for ulocation in load_init_primary_quanties()
        IMAS.IMASDD.selective_copy!(dd, dd_out, IMAS.i2p(ulocation), dd.global_time)
    end

    return dd_out
end
