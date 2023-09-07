"""
    ini_from_ods!(ini::ParametersAllInits)::IMAS.dd

The purpose of this function is to setting `ini` values based on what is in the ods
thus simplifying the logic of the init functions after it which only have to look at ini values
"""
function ini_from_ods!(ini::ParametersAllInits)::IMAS.dd

    if ini.general.init_from == :ods
        dd1 = IMAS.json2imas(ini.ods.filename)
        dd1.global_time = ini.time.simulation_start

        # equilibrium
        eqt = dd1.equilibrium.time_slice[]
        IMAS.flux_surfaces(eqt)
        if ismissing(ini.equilibrium, :R0) && !ismissing(dd1.equilibrium.vacuum_toroidal_field, :r0)
            ini.equilibrium.R0 = dd1.equilibrium.vacuum_toroidal_field.r0
        end
        if ismissing(ini.equilibrium, :B0) && !ismissing(dd1.equilibrium.vacuum_toroidal_field, :b0)
            ini.equilibrium.B0 = @ddtime dd1.equilibrium.vacuum_toroidal_field.b0
        end
        if ismissing(ini.equilibrium, :pressure_core) && !ismissing(eqt.profiles_1d, :pressure)
            ini.equilibrium.pressure_core = eqt.profiles_1d.pressure[1]
        end
        if ismissing(ini.equilibrium, :ip) && !ismissing(eqt.global_quantities, :ip)
            ini.equilibrium.ip = eqt.global_quantities.ip
        end
        if ismissing(ini.equilibrium, :xpoints)
            nx = length(eqt.boundary.x_point)
            if nx == 0
                ini.equilibrium.xpoints = :none
            elseif nx == 1
                if eqt.boundary.x_point[1].z > eqt.boundary.geometric_axis.z
                    ini.equilibrium.xpoints = :upper
                else
                    ini.equilibrium.xpoints = :lower
                end
            elseif nx == 2
                ini.equilibrium.xpoints = :double
            else
                error("cannot handle $nx x-points")
            end
        end

        # core_profiles
        if !ismissing(dd1.core_profiles.global_quantities, :ejima)
            ini.core_profiles.ejima = @ddtime(dd1.core_profiles.global_quantities.ejima)
            asdas
        end

    else
        dd1 = IMAS.dd()
    end

    return dd1
end
