#= ==================== =#
#  init equilibrium IDS  #
#= ==================== =#
"""
    init_equilibrium(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.equilibrium` starting from `ini` and `act` parameters
"""
function init_equilibrium(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    TimerOutputs.reset_timer!("init_equilibrium")
    TimerOutputs.@timeit timer "init_equilibrium" begin
        init_from = ini.general.init_from

        if init_from == :ods
            dd1 = IMAS.json2imas(ini.ods.filename)
            if !ismissing(dd1.equilibrium, :time) && length(keys(dd1.equilibrium.time)) > 0
                dd.global_time = max(dd.global_time, maximum(dd1.equilibrium.time))
                dd.equilibrium = dd1.equilibrium
                eqt = dd.equilibrium.time_slice[]
                IMAS.flux_surfaces(eqt)
            else
                init_from = :scalars
            end
        end

        boundary_from = ini.equilibrium.boundary_from

        if init_from == :scalars || (init_from == :ods && ini.equilibrium.boundary_from != :ods)

            # we make a copy because we overwrite some parameters
            # locally to this functions so that things work from
            # different entry points
            ini = deepcopy(ini)

            if init_from == :ods
                ini.equilibrium.ip = eqt.global_quantities.ip
                ini.equilibrium.R0 = dd.equilibrium.vacuum_toroidal_field.r0
                ini.equilibrium.B0 = @ddtime dd.equilibrium.vacuum_toroidal_field.b0
                ini.equilibrium.pressure_core = eqt.profiles_1d.pressure[1]

                pr, pz = eqt.boundary.outline.r, eqt.boundary.outline.z
                pr, pz = IMAS.resample_2d_path(pr, pz; n_points=101)
                pr, pz = IMAS.reorder_flux_surface!(pr, pz)
                if ini.equilibrium.boundary_from == :MXH_params
                    mxh = IMAS.MXH(pr, pz, 2)
                elseif boundary_from == :scalars
                    mxh = IMAS.MXH(
                        ini.equilibrium.R0,
                        ini.equilibrium.Z0,
                        ini.equilibrium.ϵ,
                        ini.equilibrium.κ,
                        0.0,
                        [0.0, 0.0],
                        [asin(ini.equilibrium.δ), -ini.equilibrium.ζ])
                end
                if ismissing(ini.equilibrium, :xpoints_number)
                    # if number of x-points is not set explicitly, get it from the ODS
                    ini.equilibrium.xpoints_number = length(eqt.boundary.x_point)
                end

            else
                eqt = resize!(dd.equilibrium.time_slice)
                if boundary_from == :rz_points
                    # R,Z boundary from points
                    if ismissing(ini.equilibrium, :rz_points)
                        error("ini.equilibrium.boundary_from is set as $boundary_from but rz_points wasn't set")
                    end
                    pr, pz = ini.equilibrium.rz_points[1], ini.equilibrium.rz_points[2]
                    pr, pz = IMAS.resample_2d_path(pr, pz; n_points=101)
                    pr, pz = IMAS.reorder_flux_surface!(pr, pz)
                    mxh = IMAS.MXH(pr, pz, 4)

                elseif boundary_from == :MXH_params
                    # R,Z boundary from MXH
                    if ismissing(ini.equilibrium, :MXH_params)
                        error("ini.equilibrium.boundary_from is set as $boundary_from but MXH_params wasn't set")
                    end
                    mxh = IMAS.MXH(ini.equilibrium.MXH_params)

                elseif boundary_from == :scalars
                    # R,Z boundary from scalars
                    mxh = IMAS.MXH(
                        ini.equilibrium.R0,
                        ini.equilibrium.Z0,
                        ini.equilibrium.ϵ,
                        ini.equilibrium.κ,
                        0.0,
                        [0.0, 0.0],
                        [asin(ini.equilibrium.δ), -ini.equilibrium.ζ])
                end

                # scalars consistent with MXH parametrization
                ini.equilibrium.ϵ = mxh.ϵ
                ini.equilibrium.R0 = mxh.R0
                ini.equilibrium.Z0 = mxh.Z0
                ini.equilibrium.κ = mxh.κ
                ini.equilibrium.δ = sin(mxh.s[1])
                ini.equilibrium.ζ = -mxh.s[2]
            end

            # ultimately always initialize from mxh
            init_equilibrium_boundary(eqt, mxh, ini.equilibrium.xpoints_number)

            # scalar quantities
            eqt.global_quantities.ip = ini.equilibrium.ip
            dd.equilibrium.vacuum_toroidal_field.r0 = ini.equilibrium.R0
            @ddtime dd.equilibrium.vacuum_toroidal_field.b0 = ini.equilibrium.B0

            # pressure and j_tor to be used by equilibrium solver
            if init_from == :ods
                # take p and j from ods
                eqt1 = dd1.equilibrium.time_slice[]
                eqt.profiles_1d.psi = eqt1.profiles_1d.psi
                eqt.profiles_1d.j_tor = eqt1.profiles_1d.j_tor
                eqt.profiles_1d.pressure = eqt1.profiles_1d.pressure
            else
                # guesses for pressure and j_tor
                psin = eqt.profiles_1d.psi = LinRange(0, 1, 129)
                eqt.profiles_1d.j_tor = eqt.global_quantities.ip .* (1.0 .- psin .^ 2) ./ eqt.boundary.geometric_axis.r
                eqt.profiles_1d.pressure = ini.equilibrium.pressure_core .* (1.0 .- psin)
            end

            # solve equilibrium
            act_copy = deepcopy(act)
            act_copy.ActorCHEASE.rescale_eq_to_ip = true
            ActorEquilibrium(dd, act_copy)
        end

        # field null surface
        if ini.equilibrium.field_null_surface > 0.0
            pushfirst!(dd.equilibrium.time_slice, field_null_surface(dd.equilibrium.time_slice[], ini.equilibrium.field_null_surface))
            pushfirst!(dd.equilibrium.vacuum_toroidal_field.b0, @ddtime(dd.equilibrium.vacuum_toroidal_field.b0))
            pushfirst!(dd.equilibrium.time, -Inf)
            dd.equilibrium.time_slice[1].time = -Inf
        end

        return dd
    end
end

"""
    init_equilibrium_boundary(
        eqt::IMAS.equilibrium__time_slice,
        mxh::IMAS.MXH,
        xpoints_number::Integer)

Initialize equilibrium boundary based on MXH boundary and number of x_points
"""
function init_equilibrium_boundary(
    eqt::IMAS.equilibrium__time_slice,
    mxh::IMAS.MXH,
    xpoints_number::Integer)

    if xpoints_number > 0
        # MXHboundary adds x-points
        mxhb = fitMXHboundary(mxh; lower_x_point=(xpoints_number >= 1), upper_x_point=(xpoints_number == 2))
        pr = mxhb.r_boundary
        pz = mxhb.z_boundary

        # boundary with x-points parametrized with MXH
        mxh = IMAS.MXH(pr, pz, 4)
    end
    pr, pz = mxh()

    # scalars
    eqt.boundary.minor_radius = mxh.ϵ * mxh.R0
    eqt.boundary.geometric_axis.r = mxh.R0
    eqt.boundary.geometric_axis.z = Z0 = mxh.Z0
    eqt.boundary.elongation = mxh.κ
    eqt.boundary.triangularity = sin(mxh.s[1])
    eqt.boundary.squareness = -mxh.s[2]

    # x points at maximum curvature
    x_points = []
    if xpoints_number >= 1
        i1 = argmax(abs.(IMAS.curvature(pr, pz)) .* (pz .< Z0))
        push!(x_points, (pr[i1], pz[i1]))
    end
    if xpoints_number == 2
        i2 = argmax(abs.(IMAS.curvature(pr, pz)) .* (pz .> Z0))
        push!(x_points, (pr[i2], pz[i2]))
    end
    if length(x_points) > 0
        resize!(eqt.boundary.x_point, length(x_points))
        for (k, xp_rz) in enumerate(x_points)
            eqt.boundary.x_point[k].r = xp_rz[1]
            eqt.boundary.x_point[k].z = xp_rz[2]
        end
    end

    # boundary
    eqt.boundary.outline.r = pr
    eqt.boundary.outline.z = pz

    return eqt
end

"""
    field_null_surface(eqt::IMAS.equilibrium__time_slice, scale::Real=0.5, abs_psi_boundary::Real=0.1)

Return field null surface by scaling an existing equilibrium time_slice
"""
function field_null_surface(eqt::IMAS.equilibrium__time_slice, scale::Real=0.5, abs_psi_boundary::Real=0.1)
    eqb = IMAS.equilibrium__time_slice()

    pr, pz = eqt.boundary.outline.r, eqt.boundary.outline.z
    pr, pz = IMAS.resample_2d_path(pr, pz; n_points=101)
    pr, pz = IMAS.reorder_flux_surface!(pr, pz)
    mxh = IMAS.MXH(pr, pz, 3)
    mxh.ϵ *= scale
    eqb.boundary.outline.r, eqb.boundary.outline.z = mxh()

    eqb.global_quantities.psi_boundary = abs_psi_boundary
    eqb.profiles_1d.psi = [eqb.global_quantities.psi_boundary]
    eqb.profiles_1d.f = [eqt.profiles_1d.f[end]]
    return eqb
end