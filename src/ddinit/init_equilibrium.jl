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

        # we make a copy because we overwrite some parameters
        # locally to this functions so that things work from
        # different entry points
        ini = deepcopy(ini)

        # if elongation <1.0 then expresses elongation as fraction of maximum controllable elongation estimate
        if !ismissing(ini.equilibrium, :κ) && ini.equilibrium.κ < 1.0 && !ismissing(ini.equilibrium, :ϵ)
            ini.equilibrium.κ = IMAS.elongation_limit(1.0 / ini.equilibrium.ϵ) * ini.equilibrium.κ
        end

        if init_from == :ods
            ini.equilibrium.ip = eqt.global_quantities.ip
            ini.equilibrium.R0 = dd.equilibrium.vacuum_toroidal_field.r0
            ini.equilibrium.B0 = @ddtime dd.equilibrium.vacuum_toroidal_field.b0
            ini.equilibrium.pressure_core = eqt.profiles_1d.pressure[1]

            pr, pz = eqt.boundary.outline.r, eqt.boundary.outline.z
            pr, pz = IMAS.resample_2d_path(pr, pz; n_points=101)
            pr, pz = IMAS.reorder_flux_surface!(pr, pz)

            if boundary_from == :ods
                mxh = IMAS.MXH(pr, pz, 4)
            end

            if ismissing(ini.equilibrium, :xpoints)
                # if number of x-points is not set explicitly, get it from the ODS
                n_xpoints = length(eqt.boundary.x_point)
                if n_xpoints == 0
                    ini.equilibrium.xpoints = :none
                elseif n_xpoints == 1
                    if eqt.boundary.x_point[1].z > eqt.boundary.geometric_axis.z
                        ini.equilibrium.xpoints = :upper
                    else
                        ini.equilibrium.xpoints = :lower
                    end
                elseif n_xpoints == 2
                    ini.equilibrium.xpoints = :double
                else
                    error("cannot handle $n_xpoints x-points")
                end
            end
        end

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

        # initialize dd.pulse_schedule.position_control from mxh
        init_pulse_schedule_postion_control(dd.pulse_schedule.position_control, mxh, ini.equilibrium.xpoints)

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
                cp1d.grid.psi = eqt1.profiles_1d.psi
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
        if init_from != :ods
            act_copy = deepcopy(act)
            act_copy.ActorCHEASE.rescale_eq_to_ip = true
            ActorEquilibrium(dd, act_copy)
        end

        # field null surface
        if ini.equilibrium.field_null_surface > 0.0
            field_null_surface!(dd.pulse_schedule.position_control, dd.equilibrium, ini.equilibrium.field_null_surface)
        end

        return dd
    end
end

"""
    init_pulse_schedule_postion_control(
        pc::IMAS.pulse_schedule__position_control,
        mxh::IMAS.MXH,
        xpoints::Symbol)

Initialize pulse_schedule.postion_control based on MXH boundary and number of x_points
"""
function init_pulse_schedule_postion_control(
    pc::IMAS.pulse_schedule__position_control,
    mxh::IMAS.MXH,
    xpoints::Symbol)

    if xpoints in [:lower, :upper, :double]
        # MXHboundary adds x-points
        mxhb = fitMXHboundary(mxh; lower_x_point=(xpoints in [:lower, :double]), upper_x_point=(xpoints in [:upper, :double]))
        pr = mxhb.r_boundary
        pz = mxhb.z_boundary

        # x-point information
        resize!(pc.x_point, length(mxhb.RX))
        for (k, (rx, zx)) in enumerate(zip(mxhb.RX, mxhb.ZX))
            @ddtime(pc.x_point[k].r.reference.data = rx)
            @ddtime(pc.x_point[k].r.reference.data = zx)
        end

        # boundary with x-points parametrized with MXH
        mxh = IMAS.MXH(pr, pz, 5)
    end
    pr, pz = mxh(100; adaptive=false)

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
    if xpoints in [:lower, :double]
        i1 = argmax(abs.(IMAS.curvature(pr, pz)) .* (pz .< Z0))
        push!(x_points, (pr[i1], pz[i1]))
    end
    if xpoints in [:upper, :double]
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
    for (k, (r, z)) in enumerate(zip(pr, pz))
        @ddtime(pc.boundary_outline[k].r.reference.data = r)
        @ddtime(pc.boundary_outline[k].z.reference.data = z)
    end

    return pc
end

"""
    field_null_surface!(pc::IMAS.pulse_schedule__position_control, eq::IMAS.equilibrium, scale::Real=0.5, ψp_constant::Real=0.1)

Setup field null surface as pulse_schedule.position_control.boundary_outline and insert equilibrium time slice at time=-Inf
"""
function field_null_surface!(pc::IMAS.pulse_schedule__position_control, eq::IMAS.equilibrium, scale::Real=0.75, ψp_constant::Real=0.1)
    eqt = eq.time_slice[]
    mxh = IMAS.MXH(eqt.boundary.outline.r, eqt.boundary.outline.z, 0)
    mxh.κ = 1.0
    mxh.c0 = 0.0
    mxh.ϵ *= scale
    pr, pz = mxh(length(pc.boundary_outline); adaptive=false)

    if !isempty(pc.boundary_outline) && pc.boundary_outline[1].r.reference.time[1] == -Inf
        for (k, (r, z)) in enumerate(zip(pr, pz))
            pc.boundary_outline[k].r.reference.data[1] = r
            pc.boundary_outline[k].z.reference.data[1] = z
        end
    else
        for (k, (r, z)) in enumerate(zip(pr, pz))
            pushfirst!(pc.boundary_outline[k].r.reference.time, -Inf)
            pushfirst!(pc.boundary_outline[k].r.reference.data, r)
            pushfirst!(pc.boundary_outline[k].z.reference.time, -Inf)
            pushfirst!(pc.boundary_outline[k].z.reference.data, z)
        end
    end

    if !isempty(eq.time_slice) && eq.time[1] == -Inf
        eqb = empty!(eq.time_slice[1])
    else
        eqb = IMAS.equilibrium__time_slice()
        pushfirst!(eq.time_slice, eqb)
        pushfirst!(eq.time, -Inf)
        pushfirst!(eq.vacuum_toroidal_field.b0, 0.0)
    end
    eq.vacuum_toroidal_field.b0[1] = @ddtime(eq.vacuum_toroidal_field.b0)
    eqb.time = -Inf
    eqb.global_quantities.psi_boundary = ψp_constant
    eqb.profiles_1d.psi = [ψp_constant]

    return pr, pz
end
