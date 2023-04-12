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

            # if elongation is not defined, then set it to 95% of maximum controllable elongation estimate
            if ismissing(ini.equilibrium, :κ) && !ismissing(ini.equilibrium, :ϵ)
                ini.equilibrium.κ = IMAS.elongation_limit(1.0 / ini.equilibrium.ϵ) * 0.95
            end

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

            # initialize dd.pulse_schedule.position_control from mxh
            init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxh, ini.equilibrium.xpoints_number)

            # scalar quantities
            @ddtime(dd.pulse_schedule.flux_control.i_plasma.reference.data = ini.equilibrium.ip)
            @ddtime(dd.pulse_schedule.tf.b_field_tor_vacuum_r.reference.data = ini.equilibrium.B0 * ini.equilibrium.R0)
            dd.equilibrium.vacuum_toroidal_field.r0 = ini.equilibrium.R0

            # the pressure and j_tor to be used by equilibrium solver will need to be set in dd.core_profiles
            if isempty(dd.core_profiles.profiles_1d)
                cp1d = resize!(dd.core_profiles.profiles_1d)
                if init_from == :ods
                    # take p and j from input equilibrium ods
                    eqt1 = dd1.equilibrium.time_slice[]
                    cp1d.grid.rho_tor_norm = eqt1.profiles_1d.rho_tor_norm
                    cp1d.j_tor = eqt1.profiles_1d.j_tor
                    cp1d.pressure = eqt1.profiles_1d.pressure

                else
                    # guess pressure and j_tor from input current and peak pressure
                    psin = LinRange(0, 1, 129)
                    cp1d.grid.rho_tor_norm = psin
                    cp1d.grid.psi = psin
                    cp1d.j_tor = ini.equilibrium.ip .* (1.0 .- psin .^ 2) ./ @ddtime(dd.pulse_schedule.position_control.geometric_axis.r.reference.data)
                    cp1d.pressure = ini.equilibrium.pressure_core .* (1.0 .- psin)
                end
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
    init_pulse_schedule_postion_control(
        pc::IMAS.pulse_schedule__position_control,
        mxh::IMAS.MXH,
        xpoints_number::Integer)

Initialize pulse_schedule.postion_control based on MXH boundary and number of x_points
"""
function init_pulse_schedule_postion_control(
    pc::IMAS.pulse_schedule__position_control,
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
    @ddtime(pc.minor_radius.reference.data = mxh.ϵ * mxh.R0)
    @ddtime(pc.geometric_axis.r.reference.data = mxh.R0)
    @ddtime(pc.geometric_axis.z.reference.data = mxh.Z0)
    @ddtime(pc.elongation.reference.data = mxh.κ)
    @ddtime(pc.triangularity.reference.data = sin(mxh.s[1]))
    @ddtime(pc.squareness.reference.data = -mxh.s[2])

    # x points at maximum curvature
    Z0 = mxh.Z0
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
        resize!(pc.x_point, length(x_points))
        for (k, xp_rz) in enumerate(x_points)
            @ddtime(pc.x_point[k].r.reference.data = xp_rz[1])
            @ddtime(pc.x_point[k].z.reference.data = xp_rz[2])
        end
    end

    # boundary
    resize!(pc.boundary_outline, length(pr))
    for (k,(r,z)) in enumerate(zip(pr,pz))
        @ddtime(pc.boundary_outline[k].r.reference.data = r)
        @ddtime(pc.boundary_outline[k].z.reference.data = z)
    end

    return pc
end

"""
    field_null_surface(eqt::IMAS.equilibrium__time_slice, scale::Real=0.5, abs_psi_boundary::Real=0.1)

Return field null surface by scaling an existing equilibrium time_slice
"""
function field_null_surface(eqt::IMAS.equilibrium__time_slice, scale::Real=0.5, abs_psi_boundary::Real=0.1)
    eqb = IMAS.equilibrium__time_slice()

    pr, pz = eqt.boundary.outline.r, eqt.boundary.outline.z
    pr, pz = IMAS.resample_2d_path(pr, pz; n_points=101, method=:linear)
    pr, pz = IMAS.reorder_flux_surface!(pr, pz)
    mxh = IMAS.MXH(pr, pz, 2)
    mxh.ϵ *= scale
    eqb.boundary.outline.r, eqb.boundary.outline.z = mxh()

    eqb.global_quantities.psi_boundary = abs_psi_boundary
    eqb.profiles_1d.psi = [eqb.global_quantities.psi_boundary]
    eqb.profiles_1d.f = [eqt.profiles_1d.f[end]]
    return eqb
end