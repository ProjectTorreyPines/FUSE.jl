#= ==================== =#
#  init equilibrium IDS  #
#= ==================== =#
"""
    init_equilibrium(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.equilibrium` starting from `ini` and `act` parameters
"""
function init_equilibrium(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    init_from = ini.general.init_from

    if init_from == :ods
        dd1 = IMAS.json2imas(ini.ods.filename)
        if !ismissing(dd1.equilibrium, :time) && length(keys(dd1.equilibrium.time)) > 0
            dd.global_time = max(dd.global_time, maximum(dd1.equilibrium.time))
            dd.equilibrium = dd1.equilibrium
            IMAS.flux_surfaces(dd.equilibrium.time_slice[])
        else
            init_from = :scalars
        end
    end

    if init_from == :scalars

        boundary_from = ini.equilibrium.boundary_from
        
        if boundary_from == :rz_points
            # R,Z boundary from points
            if ismissing(ini.equilibrium, :rz_points)
                error("ini.equilibrium.boundary_from is set as $boundary_from but rz_points wasn't set")
            end
            pr, pz = ini.equilibrium.rz_points[1], ini.equilibrium.rz_points[2]
            pr, pz = IMAS.resample_2d_line(pr, pz; n_points=101)
            pr, pz = IMAS.reorder_flux_surface!(pr, pz)

        else

            if boundary_from == :MXH_params
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

            # MXHboundary adds x-points
            mxhb = fitMXHboundary(mxh; lower_x_point=ini.equilibrium.x_point >= 1, upper_x_point=ini.equilibrium.x_point == 2)
            pr = mxhb.r_boundary
            pz = mxhb.z_boundary
        end

        # make arrays compatible with MXH parametrization
        pr,pz = IMAS.MXH(pr, pz, 5)()

        # ultimately always initialize from R, Z points
        init_equilibrium(
            dd.equilibrium;
            pr,
            pz,
            ini.equilibrium.x_point,
            ini.equilibrium.B0,
            ini.equilibrium.ip,
            ini.equilibrium.pressure_core)

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

"""
    init_equilibrium(
        eq::IMAS.equilibrium;
        pr::AbstractVector{<:Real},
        pz::AbstractVector{<:Real},
        x_point::Integer,
        B0::Real,
        ip::Real,
        pressure_core::Real)

Initialize equilibrium IDS based on R, Z boundary and number of x_points
"""
function init_equilibrium(
    eq::IMAS.equilibrium;
    pr::AbstractVector{<:Real},
    pz::AbstractVector{<:Real},
    x_point::Integer,
    B0::Real,
    ip::Real,
    pressure_core::Real)

    # mxh of order 2 since that's all we need to get scalars up to squareness
    mxh = IMAS.MXH(pr, pz, 2)
    R0 = mxh.R0
    Z0 = mxh.Z0
    ϵ = mxh.ϵ
    κ = mxh.κ
    δ = sin(mxh.s[1])
    ζ = -mxh.s[2]

    # scalars
    eqt = resize!(eq.time_slice)
    eqt.boundary.minor_radius = ϵ * R0
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.geometric_axis.z = Z0
    eqt.boundary.elongation = κ
    eqt.boundary.triangularity = δ
    eqt.boundary.squareness = ζ

    # x points at maximum curvature
    x_points = []
    if x_point >= 1
        i1 = argmax(abs.(IMAS.curvature(pr, pz)) .* (pz .< Z0))
        push!(x_points, (pr[i1], pz[i1]))
    end
    if x_point == 2
        i2 = argmax(abs.(IMAS.curvature(pr, pz)) .* (pz .> Z0))
        push!(x_points, (pr[i2], pz[i2]))
    end
    if length(x_points) > 0
        resize!(eqt.boundary.x_point, length(x_points))
        for (k,xp_rz) in enumerate(x_points)
            eqt.boundary.x_point[k].r = xp_rz[1]
            eqt.boundary.x_point[k].z = xp_rz[2]
        end
    end

    # scalar quantities
    eqt.global_quantities.ip = ip
    eq.vacuum_toroidal_field.r0 = R0
    @ddtime eq.vacuum_toroidal_field.b0 = B0

    # initial guesses for pressure and j_tor
    eq1d = eqt.profiles_1d
    psin = eq1d.psi = LinRange(0, 1, 129)
    eq1d.j_tor = eqt.global_quantities.ip .* (1.0 .- psin .^ 2) ./ eqt.boundary.geometric_axis.r
    eq1d.pressure = pressure_core .* (1.0 .- psin)

    # boundary
    eqt.boundary.outline.r, eqt.boundary.outline.z = pr, pz

    return eq
end

"""
    field_null_surface(eqt, scale = 0.25, abs_psi_boundary = 0.1)

Return field null surface by scaling an existing equilibrium time_slice
"""
function field_null_surface(eqt::IMAS.equilibrium__time_slice, scale::Real=0.25, abs_psi_boundary::Real=0.1)
    eqb = IMAS.equilibrium__time_slice()
    eqb.global_quantities.psi_boundary = sign(eqt.profiles_1d.psi[1] - eqt.profiles_1d.psi[end]) * abs_psi_boundary
    eqb.boundary.outline.r, eqb.boundary.outline.z, _ = IMAS.flux_surface(eqt, eqt.profiles_1d.psi[1] * (1 - scale) + eqt.profiles_1d.psi[end] * scale)
    eqb.boundary.outline.r .-= minimum(eqb.boundary.outline.r) .- minimum(IMAS.flux_surface(eqt, eqt.profiles_1d.psi[end])[1])
    eqb.profiles_1d.psi = [eqb.global_quantities.psi_boundary]
    eqb.profiles_1d.f = [eqt.profiles_1d.f[end]]
    return eqb
end