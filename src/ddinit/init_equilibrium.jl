#= ==================== =#
#  init equilibrium IDS  #
#= ==================== =#
"""
    init_equilibrium(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.equilibrium` starting from `ini` and `act` parameters
"""
function init_equilibrium(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors; set_B0_from_b_max=false, set_ip_from_q_cyl=false, kwargs...)
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
        # init equilibrium
        # use keyword argument options
        if set_B0_from_b_max
            idices = [idx for (idx,item) in enumerate(ini.build.layers) if item[1] == "hfs_TF" || item[1] == "plasma"]
            distance_tf_to_plasma = sum([item[2] for item in ini.build.layers][idices[1]+1:idices[2]-1])
            ini.equilibrium.B0 = B0_from_b_max(ini.tf.b_max, distance_tf_to_plasma, ini.equilibrium.R0 * ini.equilibrium.ϵ, ini.equilibrium.R0)
        end

        if set_ip_from_q_cyl
            ini.equilibrium.ip = ip_from_q_cyl(ini.equilibrium.R0 * ini.equilibrium.ϵ, ini.equilibrium.B0, ini.stability.q_cyl, ini.equilibrium.R0, ini.equilibrium.κ)
        end

        init_equilibrium(
            dd.equilibrium;
            B0=ini.equilibrium.B0,
            R0=ini.equilibrium.R0,
            Z0=ini.equilibrium.Z0,
            ϵ=ini.equilibrium.ϵ,
            κ=ini.equilibrium.κ,
            δ=ini.equilibrium.δ,
            ζ=ini.equilibrium.ζ,
            pressure_core=ini.equilibrium.pressure_core,
            ip=ini.equilibrium.ip,
            boundary_switch=ini.equilibrium.boundary_from,
            MXH_params=getproperty(ini.equilibrium, :MXH_params, missing),
            x_point=ini.equilibrium.x_point,
            symmetric=ini.equilibrium.symmetric)

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
        B0::Real,
        R0::Real,
        Z0::Real,
        ϵ::Real,
        κ::Real,
        δ::Real,
        ζ::Real,
        pressure_core::Real,
        ip::Real,
        boundary_switch::Symbol,
        rz_points::Union{Missing,Vector{Vector{<:Real}}}=missing,
        MXH_params::Union{Missing,Vector{<:Real}}=missing,
        x_point::Union{AbstractVector,NTuple{2},Bool}=false,
        symmetric::Bool=true)

Initialize equilibrium IDS based on some basic Miller geometry parameters
"""
function init_equilibrium(
    eq::IMAS.equilibrium;
    B0::Union{Real,Symbol},
    R0::Real,
    Z0::Real,
    ϵ::Real,
    κ::Real,
    δ::Real,
    ζ::Real,
    pressure_core::Real,
    ip::Real,
    boundary_switch::Symbol,
    rz_points::Union{Missing,Vector{Vector{<:Real}}}=missing,
    MXH_params::Union{Missing,Vector{<:Real}}=missing,
    x_point::Union{AbstractVector,NTuple{2},Bool}=false,
    symmetric::Bool=true)

    @assert -1.0 <= δ <= 1.0 "δ should be between -1.0 and 1.0"

    eqt = resize!(eq.time_slice)
    eqt.boundary.minor_radius = minor_radius = ϵ * R0
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.geometric_axis.z = Z0
    eqt.boundary.elongation = κ
    eqt.boundary.triangularity = δ
    eqt.boundary.squareness = ζ

    # x points
    if x_point === true
        mr, mz = square_miller(R0, ϵ, κ, δ, ζ; exact=true, x_points=true)
        mz .+= Z0
        i = argmax(abs.(IMAS.curvature(mr, mz)) .* (mz .< Z0))
        x_point = (mr[i], mz[i])
    end
    if isa(x_point, Union{AbstractVector,Tuple})
        resize!(eqt.boundary.x_point, 1)
        eqt.boundary.x_point[1].r = x_point[1]
        eqt.boundary.x_point[1].z = x_point[2]
        if symmetric
            resize!(eqt.boundary.x_point, 2)
            eqt.boundary.x_point[2].r = x_point[1]
            eqt.boundary.x_point[2].z = -x_point[2]
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

    if boundary_switch == :rz_points
        # R,Z boundary from: points
        if ismissing(rz_points)
            error("ini.equilibrium.boundary_from is set as $boundary_switch but rz_points wasn't set")
        end
        eqt.boundary.outline.r, eqt.boundary.outline.z = rz_points[1], rz_points[2]

    elseif boundary_switch == :MXH_params
        # R,Z boundary from: MXH
        if ismissing(MXH_params)
            error("ini.equilibrium.boundary_from is set as $boundary_switch but MXH_params wasn't set")
        end
        mxh = IMAS.MXH(MXH_params)()
        eqt.boundary.outline.r, eqt.boundary.outline.z = mxh[1], mxh[2]

    elseif boundary_switch == :scalars
        # R,Z boundary from: scalars
        eqt.boundary.outline.r, eqt.boundary.outline.z = square_miller(R0, ϵ, κ, δ, ζ; exact=true, x_points=x_point !== false)
        eqt.boundary.outline.z .+= Z0
    end

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