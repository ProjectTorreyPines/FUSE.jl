import SpecialFunctions
import QuadGK
import Interpolations
import PolygonOps
import LibGEOS
import GeoInterface

#= =============== =#
#  Shape functions  #
#= =============== =#

function layer_shape_message(shape_function_index)
    return "layer.shape=$(shape_function_index) is invalid. Valid options are:
            -2: Offset & convex-hull
            -1: Offset
             1: Princeton D exact  (shape_parameters = [])
             2: Princeton D approx (shape_parameters = [])
             3: Princeton D scaled (shape_parameters = [height])
             4: rectangle          (shape_parameters = [height])
             5: double_ellipse     (shape_parameters = [centerpost_height, height])
             6: rectangle_ellipse  (shape_parameters = [height])
             7: tripple-arc        (shape_parameters = [height, small_radius, mid_radius, small_coverage, mid_coverage])
             8: miller             (shape_parameters = [elongation, triangularity])
             9: square_miller      (shape_parameters = [elongation, triangularity, squareness])
            10: spline             (shape_parameters = [hfact, rz...)
            11: racetrack          (shape_parameters = [height, radius_ratio])
            12: silo               (shape_parameters = [h_start, h_end)
           10x: shape + z_offset   (shape_parameters = [..., z_offset])
          100x: negative shape     (shape_parameters = [...])"
end

function initialize_shape_parameters(shape_function_index, r_obstruction, z_obstruction, r_start, r_end, clearance)
    height = maximum(z_obstruction) - minimum(z_obstruction) + clearance * 2.0
    z_offset = (maximum(z_obstruction) + minimum(z_obstruction)) / 2.0
    shape_parameters = nothing
    if shape_function_index ∈ (Int(_offset_), Int(_negative_offset_), Int(_convex_hull_))
        return nothing
    else
        shape_index_mod = shape_function_index
        is_negative_D = false
        if shape_index_mod > 1000
            shape_index_mod = mod(shape_function_index, 1000)
            is_negative_D = true
        end
        is_z_offset = false
        if shape_index_mod > 100
            shape_index_mod = mod(shape_function_index, 100)
            is_z_offset = true
        end
        if shape_index_mod == Int(_princeton_D_)
            shape_parameters = Float64[]
        elseif shape_index_mod == Int(_princeton_D_scaled_)
            shape_parameters = [height]
        elseif shape_index_mod == Int(_double_ellipse_)
            r_center = (maximum(r_obstruction) + minimum(r_obstruction)) / 2.0
            centerpost_height = (maximum(z_obstruction) - minimum(z_obstruction)) * 2.0 / 3.0
            shape_parameters = [r_center, centerpost_height, height]
        elseif shape_index_mod == Int(_circle_ellipse_)
            centerpost_height = (maximum(z_obstruction) - minimum(z_obstruction)) * 2.0 / 3.0
            shape_parameters = [centerpost_height, height]
        elseif shape_index_mod == Int(_rectangle_)
            shape_parameters = [height]
        elseif shape_index_mod == Int(_racetrack_)
            shape_parameters = [height, 0.25]
        elseif shape_index_mod == Int(_triple_arc_)
            shape_parameters = [height, 0.5, 0.5, 45, 45]
        elseif shape_index_mod ∈ (Int(_miller_), Int(_square_miller_))
            mxh = IMAS.MXH(r_obstruction, z_obstruction, 2)
            shape_parameters = [mxh.κ, sin(mxh.s[1])]
            if shape_index_mod == Int(_square_miller_)
                push!(shape_parameters, -mxh.s[2])
            end
        elseif shape_index_mod == Int(_spline_)
            n = 1
            R = range(r_start, r_end, 2 + n)[2:end-1]
            Z = range(height / 2.0, height / 2.0, 2 + n)[2:end-1]
            shape_parameters = Float64[0.8]
            for (r, z) in zip(R, Z)
                append!(shape_parameters, [r, z])
            end
        elseif shape_index_mod == Int(_silo_)
            shape_parameters = [height, height / 2.0]
        end
    end
    if shape_parameters === nothing
        error(shape_function_index)
    end
    if is_z_offset
        push!(shape_parameters, z_offset)
    end
    return shape_parameters
end

function shape_function(shape_function_index::Int; resolution::Float64)
    func = nothing
    if shape_function_index ∈ (Int(_offset_), Int(_negative_offset_), Int(_convex_hull_))
        return nothing
    else
        shape_index_mod = shape_function_index
        is_negative_D = false
        if shape_index_mod > 1000
            shape_index_mod = mod(shape_function_index, 1000)
            is_negative_D = true
        end
        is_z_offset = false
        if shape_index_mod > 100
            shape_index_mod = mod(shape_function_index, 100)
            is_z_offset = true
        end
        if shape_index_mod == Int(_princeton_D_)
            func = princeton_D_approx
        elseif shape_index_mod == Int(_princeton_D_scaled_)
            func = princeton_D_scaled
        elseif shape_index_mod == Int(_circle_ellipse_)
            func = circle_ellipse
        elseif shape_index_mod == Int(_double_ellipse_)
            func = double_ellipse
        elseif shape_index_mod == Int(_rectangle_)
            func = rectangle_shape
        elseif shape_index_mod == Int(_racetrack_)
            func = racetrack
        elseif shape_index_mod == Int(_triple_arc_)
            func = triple_arc
        elseif shape_index_mod == Int(_miller_)
            func = miller_Rstart_Rend
        elseif shape_index_mod == Int(_square_miller_)
            func = square_miller_Rstart_Rend
        elseif shape_index_mod == Int(_spline_)
            func = spline_shape
        elseif shape_index_mod == Int(_silo_)
            func = silo
        end
    end
    if func === nothing
        error(layer_shape_message(shape_index_mod))
    end

    # resolution
    rfunc(args...) = func(args...; resolution)

    # zoffset
    zfunc = rfunc
    if is_z_offset
        zfunc(args...) = begin
            R, Z = rfunc(args[1:end-1]...)
            Z .+= args[end]
            return R, Z
        end
    end

    # neg-D
    dfunc = zfunc
    if is_negative_D
        dfunc(args...) = begin
            R, Z = zfunc(args...)
            R = -(R .- args[1]) .+ args[2]
            return reverse(R), reverse(Z)
        end
    end

    # uniform resampling
    myfunc = dfunc
    function resampled_myfunc(args...; resample=true)
        if resample
            return IMAS.resample_2d_path(myfunc(args...)...; method=:linear)
        else
            return myfunc(args...)
        end
    end

    return resampled_myfunc
end

"""
    optimize_outline(r_obstruction, z_obstruction, hfs_thickness, lfs_thickness, func, r_start, r_end, shape_parameters; verbose=false)

Find shape parameters that generate smallest shape and target clearance from an obstruction
"""
function optimize_outline(
    r_obstruction::Vector{Float64},
    z_obstruction::Vector{Float64},
    hfs_thickness::Float64,
    lfs_thickness::Float64,
    func::Function,
    r_start::Float64,
    r_end::Float64,
    shape_parameters::Vector{Float64};
    use_curvature::Bool=true,
    verbose::Bool=false)

    rz_obstruction = collect(zip(r_obstruction, z_obstruction))
    initial_guess = deepcopy(shape_parameters)

    if length(shape_parameters) in (0, 1)
        func(r_start, r_end, shape_parameters...)

    else
        function cost_shape(
            r_obstruction::Vector{Float64},
            z_obstruction::Vector{Float64},
            rz_obstruction::Vector{Tuple{Float64,Float64}},
            hfs_thickness::Float64,
            lfs_thickness::Float64,
            func::Function,
            r_start::Float64,
            r_end::Float64,
            shape_parameters::Vector{Float64};
            use_curvature::Bool,
            verbose::Bool=false)

            #@show r_start, r_end, shape_parameters

            R0, Z0 = func(r_start, r_end, shape_parameters...)

            if hfs_thickness == 0 || lfs_thickness == 0
                target_clearance = (hfs_thickness + lfs_thickness) / 2.0
                R, Z = R0, Z0
                hbuf = 0.0
                lbuf = 0.0
            else
                target_clearance = min(hfs_thickness, lfs_thickness)
                if hfs_thickness != lfs_thickness
                    hbuf = hfs_thickness - target_clearance
                    lbuf = lfs_thickness - target_clearance
                    R, Z = buffer(R0, Z0, -hbuf, -lbuf)
                else
                    R, Z = R0, Z0
                    hbuf = 0.0
                    lbuf = 0.0
                end
            end

            # R1, Z1 = buffer(R, Z, -target_clearance)
            # plot()
            # plot!(R0,Z0,label="R0,Z0")
            # plot!(R,Z;aspect_ratio=:equal,label="R,Z")
            # vline!([r_start+hbuf, r_end-lbuf];primary=false)
            # plot!(R1, Z1,label="R1,Z1",ls=:dash)
            # plot!()
            # display(plot!(r_obstruction,z_obstruction;color=:black,label="obstruction"))

            # disregard near r_start and r_end where optimizer has no control and shape is allowed to go over obstruction
            index = (.!)(isapprox.(R, r_start + hbuf + target_clearance) .|| isapprox.(R, r_end - lbuf - target_clearance))
            Rv = view(R, index)
            Zv = view(Z, index)

            # no polygon crossings  O(N)
            cost_inside = 0.0
            for (r, z) in zip(Rv, Zv)
                inpoly = PolygonOps.inpolygon((r, z), rz_obstruction)
                cost_inside += inpoly
            end

            # target clearance  O(1)
            minimum_distance, mean_distance = IMAS.min_mean_distance_two_shapes(Rv, Zv, r_obstruction, z_obstruction)
            cost_mean_distance = (mean_distance - target_clearance) / ((hfs_thickness + lfs_thickness) / 2.0)
            if minimum_distance < target_clearance
                cost_min_clearance = (minimum_distance - target_clearance) / target_clearance
            else
                cost_min_clearance = 0.0
            end

            # curvature
            cost_max_curvature = 0.0
            if use_curvature
                curvature = abs.(IMAS.curvature(R, Z))
                cost_max_curvature = atan(1.0 / (1.0 - maximum(curvature)^2)) / pi * 2
            end

            if verbose
                @show target_clearance
                @show minimum_distance
                @show cost_min_clearance^2
                @show cost_mean_distance^2
                @show cost_inside^2
                @show cost_max_curvature^2
                println()
            end

            # return cost
            return norm((cost_min_clearance, cost_mean_distance, cost_inside, cost_max_curvature))
        end

        initial_guess = copy(shape_parameters)
        res = Optim.optimize(
            shape_parameters -> cost_shape(r_obstruction, z_obstruction, rz_obstruction, hfs_thickness, lfs_thickness, func, r_start, r_end, shape_parameters; use_curvature),
            initial_guess, length(shape_parameters) == 1 ? Optim.BFGS() : Optim.NelderMead(), Optim.Options(; iterations=10000, f_tol=1E-4, x_tol=1E-3))
        shape_parameters = Optim.minimizer(res)
        if verbose
            cost_shape(r_obstruction, z_obstruction, rz_obstruction, hfs_thickness, lfs_thickness, func, r_start, r_end, shape_parameters; use_curvature, verbose=true)
            println(res)
        end
    end

    # R, Z = func(r_start, r_end, shape_parameters...; resample=false)
    # plot(func(r_start, r_end, initial_guess...); markershape=:x, label="initial guess")
    # plot!(r_obstruction, z_obstruction, ; markershape=:x, label="obstruction")
    # display(plot!(R, Z; markershape=:x, aspect_ratio=:equal, label="final"))

    return shape_parameters
end

"""
    struveL(ν, z)

Modified Struve function.
Ported from https://github.com/gwater/Struve.jl
(Struve had some requirements conflicts with FUSE)
"""
function struveL(ν, z)
    integrate(f, a, b) = QuadGK.quadgk(f, a, b)[1]
    _M_integral(ν, z) = -2(0.5z)^ν / (sqrt(pi) * SpecialFunctions.gamma(ν + 0.5)) * integrate(t -> exp(-z * t) * (1 - t^2)^(ν - 0.5), 0, 1)
    return SpecialFunctions.besseli(ν, z) + _M_integral(ν, z)
end

"""
    ellipse(a::T, b::T, t0::T, t1::T, x0::T, z0::T; n_points::Integer=100) where {T<:Real}

Simple ellipse shape function
"""
function ellipse(a::T, b::T, t0::T, t1::T, x0::T, z0::T; n_points::Integer=100) where {T<:Real}
    t = range(t0, t1, n_points)
    x = a .* cos.(t) .+ x0
    z = b .* sin.(t) .+ z0
    return x, z
end

"""
    princeton_D_approx(r_start::T, r_end::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

An pproximate version of the "Princeton-D" constant tension shape for a TF coil that is built with ellipses
References: Gralnick, S. L.; Tenney, F. H. Analytic Solutions for Constant‐tension Coil Shapes. J. Appl. Phys. 1976, 47, 7
"""
function princeton_D_approx(r_start::T, r_end::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    n_points = Int(round(n_points * resolution))

    r1 = r_start
    r2 = r_end
    k = 0.5 * log(r2 / r1)
    r0 = sqrt(r1 * r2)

    # analytic equations for coil parameters
    centerpost_maxz = 2 * pi * r0 * k * SpecialFunctions.besseli(1, k) / 2 # Gralnick Eq. 28
    coil_maxz = pi * r0 * k * (SpecialFunctions.besseli(1, k) + struveL(1, k) + 2 / pi) / 2  # Gralnick Eq. 34

    # number of points
    ra_e1 = r0 - r1
    zb_e1 = coil_maxz - centerpost_maxz
    ra_e2 = r2 - r0
    zb_e2 = coil_maxz
    e1_e2_ratio = sqrt(ra_e1^2 + zb_e1^2) / sqrt(ra_e2^2 + zb_e2^2)
    n1_points = max(3, Int(round(n_points * resolution * e1_e2_ratio)))
    n2_points = Int(round(n_points * resolution))

    # make ellipses to approximate equal-tension arc
    r1, z1 = ellipse(r0 - r1, coil_maxz - centerpost_maxz, float(π), float(π / 2), r0, centerpost_maxz; n_points=n1_points)
    r2, z2 = ellipse(r2 - r0, coil_maxz, float(π / 2), float(0.0), r0, 0.0; n_points=n2_points)

    # stack
    R = [r1[1:end-1]; r2[1:end-1]; r2[end:-1:2]; r1[end:-1:1]; r1[1]]
    Z = [z1[1:end-1]; z2[1:end-1]; -z2[end:-1:2]; -z1[end:-1:1]; z1[1]]

    return R, Z
end

"""
    princeton_D_scaled(r_start::T, r_end::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

This routine calculates a "shortened" TF coil shape that foregoes the equal-tension "Princeton-Dee" for
a squater, more space-efficient shape. It replicates the inboard curve of the equal-tension arc, but
decreases the height of the coil to match a given value.
"""
function princeton_D_scaled(r_start::T, r_end::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    n_points = Int(round(n_points * resolution))

    r1 = r_start
    r2 = r_end
    k = 0.5 * log(r2 / r1)
    r0 = sqrt(r1 * r2)

    centerpost_maxz = 2 * pi * r0 * k * SpecialFunctions.besseli(1, k) / 2 # Gralnick Eq. 28
    coil_maxz = pi * r0 * k * (SpecialFunctions.besseli(1, k) + struveL(1, k) + 2 / pi) / 2  # Gralnick Eq. 34

    centerpost_height = height - (coil_maxz - centerpost_maxz) * 2

    inboard_curve_dz = coil_maxz - centerpost_maxz
    centerpost_maxz = centerpost_height / 2
    coil_maxz = centerpost_maxz + inboard_curve_dz

    # number of points
    ra_e1 = r0 - r1
    zb_e1 = coil_maxz - centerpost_maxz
    ra_e2 = r2 - r0
    zb_e2 = coil_maxz
    e1_e2_ratio = sqrt(ra_e1^2 + zb_e1^2) / sqrt(ra_e2^2 + zb_e2^2)
    n1_points = max(3, Int(round(n_points * resolution * e1_e2_ratio)))
    n2_points = Int(round(n_points * resolution))

    # make ellipses to approximate equal-tension arc
    r1, z1 = ellipse(r0 - r1, coil_maxz - centerpost_maxz, float(π), float(π / 2), r0, centerpost_maxz; n_points=n1_points)
    r2, z2 = ellipse(r2 - r0, coil_maxz, float(π / 2), float(0.0), r0, 0.0; n_points=n2_points)

    # stack
    R = [r1[1:end-1]; r2[1:end-1]; r2[end:-1:2]; r1[end:-1:1]; r1[1]]
    Z = [z1[1:end-1]; z2[1:end-1]; -z2[end:-1:2]; -z1[end:-1:1]; z1[1]]

    return R, Z
end

"""
    double_ellipse(r_start::T, r_end::T, r_center::T, centerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

double ellipse shape
"""
function double_ellipse(r_start::T, r_end::T, r_center::T, centerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    height = abs(height)
    centerpost_height = mirror_bound(abs(centerpost_height), 0.0, height)
    r_center = mirror_bound(r_center, r_start, r_end)
    return double_ellipse(r_start, r_end, r_center, centerpost_height, 0.0, height; n_points, resolution)
end

function double_ellipse(r_start::T, r_end::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    return double_ellipse(r_start, r_end, (r_start + r_end) / 2, height, 0.0, height; n_points, resolution)
end

function double_ellipse(r_start::T, r_end::T, r_center::T, centerpost_height::T, outerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    # inner
    r_e1 = r_center
    z_e1 = centerpost_height / 2.0
    ra_e1 = r_center - r_start
    zb_e1 = (height - centerpost_height) / 2.0

    # outer
    r_e2 = r_center
    z_e2 = outerpost_height / 2.0
    ra_e2 = r_end - r_center
    zb_e2 = (height - outerpost_height) / 2.0

    # number of points
    e1_e2_ratio = sqrt(ra_e1^2 + zb_e1^2) / sqrt(ra_e2^2 + zb_e2^2)
    n1_points = max(3, Int(round(n_points * resolution * e1_e2_ratio)))
    n2_points = Int(round(n_points * resolution))

    # ellipses
    r1, z1 = ellipse(ra_e1, zb_e1, float(π), float(π / 2), r_e1, z_e1; n_points=n1_points)
    r2, z2 = ellipse(ra_e2, zb_e2, float(π / 2), float(0.0), r_e2, z_e2; n_points=n2_points)

    # stack
    R = [r1[1:end-1]; r2[1:end-1]; r2[end:-1:2]; r1[end:-1:1]; r1[1]]
    Z = [z1[1:end-1]; z2[1:end-1]; -z2[end:-1:2]; -z1[end:-1:1]; z1[1]]

    return R, Z
end

"""
    circle_ellipse(r_start::T, r_end::T, centerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

circle ellipse shape (parametrization of TF coils used in GATM)

Special case of the double ellipse shape, where the inner ellipse is actually a circle
"""
function circle_ellipse(r_start::T, r_end::T, centerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    height = abs(height)
    centerpost_height = mirror_bound(abs(centerpost_height), 0.0, height)
    r_circle = (height - centerpost_height) / 2.0
    r_center = r_start + r_circle
    return double_ellipse(r_start, r_end, r_center, centerpost_height, 0.0, height; n_points, resolution)
end

"""
    rectangle_shape(r_start::T, r_end::T, z_low::T, z_high::T; n_points::Int=400, resolution::Float64=0.0) where {T<:Real}

Asymmetric rectangular shape

NOTE: by default resolution=0, which returns 5 points
"""
function rectangle_shape(r_start::T, r_end::T, z_low::T, z_high::T; n_points::Int=400, resolution::Float64=0.0) where {T<:Real}
    if resolution == 0.0
        R = [r_start, r_start, r_end, r_end, r_start]
        Z = [z_low, z_high, z_high, z_low, z_low]
    else
        n_points = Int(round(n_points / 4 * resolution))
        R = vcat(
            range(r_start, r_start, n_points),
            range(r_start, r_end, n_points)[2:end],
            range(r_end, r_end, n_points)[2:end],
            range(r_end, r_start, n_points)[2:end])
        Z = vcat(
            range(z_low, z_high, n_points),
            range(z_high, z_high, n_points)[2:end],
            range(z_high, z_low, n_points)[2:end],
            range(z_low, z_low, n_points)[2:end])
    end
    return R, Z
end

"""
    rectangle_shape(r_start::T, r_end::T, height::T; n_points::Integer=400, resolution::Float64=0.0) where {T<:Real}

Symmetric rectangular shape

NOTE: by default resolution=0, which returns 5 points
"""
function rectangle_shape(r_start::T, r_end::T, height::T; n_points::Integer=400, resolution::Float64=0.0) where {T<:Real}
    Δ = height / 2.0
    return rectangle_shape(r_start, r_end, -Δ, Δ; n_points, resolution)
end

"""
    racetrack(r_start::T, r_end::T, z_low::T, z_high::T, radius_ratio::T; n_points::Int=400, resolution::Float64=1.0) where {T<:Real}

Asymmetric racetrack shape
"""
function racetrack(r_start::T, r_end::T, z_low::T, z_high::T, radius_ratio::T; n_points::Int=400, resolution::Float64=1.0) where {T<:Real}
    R, Z = rectangle_shape(r_start, r_end, z_low, z_high; n_points, resolution)
    if radius_ratio != 0.0
        radius_ratio = mirror_bound(abs(radius_ratio), 0.0, 1.0 - 1E-3)
        radius = min((z_high - z_low), (r_end - r_start)) / 2.0 * radius_ratio
        R, Z = buffer(R, Z, -radius)
        R, Z = buffer(R, Z, +radius)
    end
    return R, Z
end

"""
    racetrack(r_start::T, r_end::T, height::T, radius_ratio::T; n_points::Integer=400, resolution::Float64=1.0) where {T<:Real}

Symmetric racetrack shape
"""
function racetrack(r_start::T, r_end::T, height::T, radius_ratio::T; n_points::Integer=400, resolution::Float64=1.0) where {T<:Real}
    Δ = height / 2.0
    return racetrack(r_start, r_end, -Δ, Δ, radius_ratio; n_points, resolution)
end

"""
    triple_arc(
        r_start::T,
        r_end::T,
        height::T,
        small_radius::T,
        mid_radius::T,
        small_coverage::T,
        mid_coverage::T;
        min_small_radius_fraction::T=0.2,
        min_mid_radius_fraction::T=min_small_radius_fraction * 2.0,
        n_points::Integer=100,
        resolution::Float64=1.0) where {T<:Real}

TrippleArc shape. Angles are in degrees.
"""
function triple_arc(
    r_start::T,
    r_end::T,
    height::T,
    small_radius::T,
    mid_radius::T,
    small_coverage::T,
    mid_coverage::T;
    n_points::Integer=400,
    resolution::Float64=1.0) where {T<:Real}

    n_points = Int(round(n_points / 4.0 * resolution))

    height = abs(height) / 2.0
    small_radius = mirror_bound(small_radius, 0.01, 1.0) * height
    mid_radius = mirror_bound(mid_radius, small_radius / height, 1.0) * height
    small_coverage = mirror_bound(small_coverage, 0.0, 180.0) * pi / 180
    mid_coverage = mirror_bound(mid_coverage, 0.0, 180.0) * pi / 180

    asum = small_coverage + mid_coverage

    # small arc
    theta = range(0, small_coverage, n_points)
    small_arc_R = r_start .+ small_radius .* (1 .- cos.(theta))
    small_arc_Z = height .+ small_radius .* sin.(theta)

    # mid arc
    theta = range(small_coverage, asum, n_points)
    mid_arc_R = small_arc_R[end] .+ mid_radius .* (cos.(small_coverage) .- cos.(theta))
    mid_arc_Z = small_arc_Z[end] .+ mid_radius .* (sin.(theta) .- sin.(small_coverage))

    # large arc
    theta = range(theta[end], pi, n_points)
    large_radius = mid_arc_Z[end] / sin(pi - asum)
    large_arc_R = mid_arc_R[end] .+ large_radius .* (cos.(pi .- theta) .- cos.(pi .- asum))
    large_arc_Z = mid_arc_Z[end] .- large_radius .* (sin(asum) .- sin.(pi .- theta))

    R = vcat(small_arc_R, mid_arc_R[2:end], large_arc_R[2:end])
    R = vcat(R, reverse(R)[2:end])
    Z = vcat(small_arc_Z, mid_arc_Z[2:end], large_arc_Z[2:end])
    Z = vcat(Z, -reverse(Z)[2:end])

    # Add vertical
    R = vcat(range(r_start, r_start, n_points), R)
    Z = vcat(range(-height, height, n_points), Z)

    # Resize to ensure r_start to r_end
    factor = (r_end - r_start) / (maximum(R) - minimum(R))
    Z = Z .* factor
    R = (R .- minimum(R)) .* factor .+ r_start

    return R, Z
end

"""
    miller(R0::T, rmin_over_R0::T, elongation::T, triangularity::T; n_points::Integer=201, resolution::Float64=1.0) where {T<:Real}

Miller shape
"""
function miller(R0::T, rmin_over_R0::T, elongation::T, triangularity::T; n_points::Integer=201, resolution::Float64=1.0) where {T<:Real}
    n_points = Int(floor(n_points * resolution / 2)) * 2 + 1

    θ = range(0, 2π, n_points)
    triangularity = mirror_bound(triangularity, -1.0, 1.0)
    δ₀ = asin(triangularity)
    R = R0 * (1 .+ rmin_over_R0 .* cos.(θ .+ δ₀ * sin.(θ)))
    Z = R0 * (rmin_over_R0 * elongation * sin.(θ))
    R[end] = R[1]
    Z[end] = Z[1]
    return R, Z
end

"""
    miller_Rstart_Rend(r_start::T, r_end::T, elongation::T, triangularity::T; n_points::Int=201, resolution::Float64=1.0) where {T<:Real}

Miller shape
"""
function miller_Rstart_Rend(r_start::T, r_end::T, elongation::T, triangularity::T; n_points::Int=201, resolution::Float64=1.0) where {T<:Real}
    elongation = mirror_bound(abs(elongation), 1.0, 3.0)
    triangularity = mirror_bound(triangularity, -1.0, 1.0)
    return miller((r_end + r_start) / 2.0, (r_end - r_start) / (r_end + r_start), elongation, triangularity; n_points, resolution)
end

"""
    square_miller(
        R0::T,
        rmin_over_R0::T,
        elongation::T,
        triangularity::T,
        squareness::T;
        upper_x_point::Bool=false,
        lower_x_point::Bool=false,
        exact::Bool=false,
        n_points::Integer=201,
        resolution::Float64=1.0) where {T<:Real}

Miller contour with squareness (via MXH parametrization)

`exact=true` optimizes diverted shape to match desired Miller elongation, triangularity, squareness
"""
function square_miller(
    R0::T,
    rmin_over_R0::T,
    elongation::T,
    triangularity::T,
    squareness::T;
    upper_x_point::Bool=false,
    lower_x_point::Bool=false,
    exact::Bool=false,
    n_points::Integer=201,
    resolution::Float64=1.0) where {T<:Real}

    n_points = Int(floor(n_points * resolution / 2)) * 2 + 1

    mxh = IMAS.MXH(R0, 2)
    mxh.ϵ = rmin_over_R0
    mxh.κ = elongation
    mxh.s[1] = asin(triangularity)
    mxh.s[2] = -squareness

    if exact
        func = fitMXHboundary
    else
        func = MXHboundary
    end
    mxhb = func(mxh; upper_x_point, lower_x_point, n_points)

    return mxhb.r_boundary, mxhb.z_boundary
end

"""
    square_miller_Rstart_Rend(r_start::Real, r_end::Real, elongation::Real, triangularity::Real, squareness::Real; n_points::Int=401, resolution::Float64=1.0)

Miller with squareness contour
"""
function square_miller_Rstart_Rend(r_start::Real, r_end::Real, elongation::Real, triangularity::Real, squareness::Real; n_points::Int=401, resolution::Float64=1.0)
    elongation = mirror_bound(abs(elongation), 1.0, 3.0)
    triangularity = mirror_bound(triangularity, -1.0, 1.0)
    squareness = mirror_bound(abs(squareness), 0.0, 1.0)
    return square_miller((r_end + r_start) / 2.0, (r_end - r_start) / (r_end + r_start), elongation, triangularity, squareness; n_points, resolution)
end

function spline_shape(r::Vector{T}, z::Vector{T}; n_points::Int=101, resolution::Float64=1.0) where {T<:Real}
    n_points = Int(round(n_points * resolution))

    r = vcat(r[1], r[1], r, r[end], r[end])
    z = vcat(0, z[1] / 2, z, z[end] / 2, 0)
    d = cumsum(sqrt.(vcat(0, diff(r)) .^ 2.0 .+ vcat(0, diff(z)) .^ 2.0))

    itp_r = Interpolations.interpolate(d, r, Interpolations.FritschButlandMonotonicInterpolation())
    itp_z = Interpolations.interpolate(d, z, Interpolations.FritschButlandMonotonicInterpolation())

    D = range(d[1], d[end], n_points)
    R, Z = itp_r.(D), itp_z.(D)
    R[end] = R[1]
    Z[end] = Z[1]
    return R, Z
end

"""
    spline_shape(r_start::T, r_end::T, hfact::T, rz...; n_points::Integer=101, resolution::Float64=1.0) where {T<:Real}

Spline shape
"""
function spline_shape(r_start::T, r_end::T, hfact::T, rz...; n_points::Integer=101, resolution::Float64=1.0) where {T<:Real}
    rz = collect(rz)
    R = rz[1:2:end]
    Z = rz[2:2:end]
    hfact_max = 0.5 + (minimum(R) - r_start) / (r_end - r_start) / 2.0
    hfact = min(abs(hfact), hfact_max)
    h = maximum(Z) * hfact
    r = vcat(r_start, R, r_end, reverse(R), r_start)
    z = vcat(h, Z, 0, -reverse(Z), -h)
    return spline_shape(r, z; n_points, resolution)
end

"""
    xy_polygon(x::T, y::T) where {T<:AbstractVector{<:Real}}

Returns LibGEOS.Polygon from x and y arrays
"""
function xy_polygon(x::T, y::T) where {T<:AbstractVector{<:Real}}
    if (x[1] ≈ x[end]) && (y[1] ≈ y[end])
        coords = [[[x[i], y[i]] for i in 1:length(x)]]
        coords[1][end] .= coords[1][1]
    else
        coords = [[i > length(x) ? [x[1], y[1]] : [x[i], y[i]] for i in 1:length(x)+1]]
    end

    if (coords[1][1][1] != coords[1][end][1]) && (coords[1][1][1] ≈ coords[1][end][1])
        coords[1][end][1] = coords[1][1][1]
    end
    if (coords[1][1][2] != coords[1][end][2]) && (coords[1][1][2] ≈ coords[1][end][2])
        coords[1][end][2] = coords[1][1][2]
    end
    return LibGEOS.Polygon(coords)
end

"""
    xy_polygon(coords::T) where {T<:Vector{Tuple{Float64, Float64}}}

Returns LibGEOS.Polygon from vector of tuple with x, y points
"""
function xy_polygon(coords::T) where {T<:Vector{Tuple{Float64,Float64}}}
    return LibGEOS.Polygon([[[coord[1], coord[2]] for coord in coords]])
end

"""
    xy_polygon(layer::Union{IMAS.build__layer,IMAS.build__structure})

Returns LibGEOS.Polygon from IMAS build layer or structure
"""
function xy_polygon(layer::Union{IMAS.build__layer,IMAS.build__structure})
    return xy_polygon(layer.outline.r, layer.outline.z)
end

"""
    buffer(x::AbstractVector{T}, y::AbstractVector{T}, b::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}

Buffer polygon defined by x,y arrays by a quantity b
"""
function buffer(x::AbstractVector{T}, y::AbstractVector{T}, b::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}
    poly = xy_polygon(x, y)
    poly_b = LibGEOS.buffer(poly, b)
    coords = GeoInterface.coordinates(poly_b)[1]
    x_b = T[v[1] for v in coords]
    y_b = T[v[2] for v in coords]
    return x_b, y_b
end

"""
    buffer(x::AbstractVector{T}, y::AbstractVector{T}, b_hfs::T, b_lfs::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}

Buffer polygon defined by x,y arrays by a quantity b_hfs to the left and b_lfs to the right
"""
function buffer(x::AbstractVector{T}, y::AbstractVector{T}, b_hfs::T, b_lfs::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}
    x_b, y_b = buffer(x, y, 0.5 * (b_lfs + b_hfs))
    x_offset = 0.5 * (b_lfs - b_hfs)
    x_b .+= x_offset
    return x_b, y_b
end

"""
    limit_curvature(x::AbstractVector{T}, y::AbstractVector{T}, max_curvature::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}

Limit maximum curvature of a polygon described by x,y arrays
"""
function limit_curvature(x::AbstractVector{T}, y::AbstractVector{T}, max_curvature::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}
    @assert max_curvature > 0.0
    poly = xy_polygon(x, y)
    poly_b = LibGEOS.buffer(LibGEOS.buffer(poly, -max_curvature), max_curvature)
    coords = GeoInterface.coordinates(poly_b)[1]
    x_b = T[v[1] for v in coords]
    y_b = T[v[2] for v in coords]
    return x_b, y_b
end

"""
    volume_no_structures(layer::IMAS.build__layer, structures::IMAS.IDSvector{<:IMAS.build__structure})

Returns volume of the layer without structures
"""
function volume_no_structures(layer::IMAS.build__layer, structures::IMAS.IDSvector{<:IMAS.build__structure})
    vol = 0.0
    for structure in structures
        vol += layer_structure_intersect_volume(layer, structure)
    end
    return layer.volume - vol
end

# this expressions is added here because volume_no_structures is not a IMAS function
IMAS.dynamic_expressions["build.layer[:].volume_no_structures"] =
    (; build, layer, _...) -> volume_no_structures(layer, build.structure)

"""
    layer_structure_intersect_volume(layer::IMAS.build__layer, structure::IMAS.build__structure)

Returns volume of the intersection between build layer volume and structure volume
"""
function layer_structure_intersect_volume(layer::IMAS.build__layer, structure::IMAS.build__structure)
    if layer.type ∈ (Int(_in_), Int(_out_))
        return layer.volume
    elseif layer.type ∈ (Int(_tf_), Int(_plasma_))
        return layer.volume
    elseif layer.side ∈ (Int(_hfs_), Int(_lfs_))
        i = IMAS.index(layer)
        if layer.side == Int(_hfs_)
            layer_in = IMAS.parent(layer)[i+1]
        else
            layer_in = IMAS.parent(layer)[i-1]
        end
    end
    layer_poly = xy_polygon(layer)
    layer_in_poly = xy_polygon(layer_in)
    ring_poly = LibGEOS.difference(layer_poly, layer_in_poly)
    structure_poly = xy_polygon(structure)

    if structure.toroidal_extent == 2 * pi
        toroidal_angles = [0.0]
    else
        toroidal_angles = structure.toroidal_angles
    end

    vol = 0.0
    for poly in GeoInterface.coordinates(LibGEOS.intersection(ring_poly, structure_poly))
        pr = [v[1] for v in poly]
        pz = [v[2] for v in poly]
        vol += IMAS.area(pr, pz) * structure.toroidal_extent * length(toroidal_angles)
    end

    return vol
end

"""
    P_LH_threshold_from_scalars(Bt0::Real, nel::Real,surface_area::Real)

Calculates the L to H transition threshold according to 2008 scaling law
"""
function P_LH_threshold_from_scalars(Bt0::Real, nel::Real, surface_area::Real)
    return 0.049 * abs(Bt0)^0.8 * nel^0.72 * surface_area^0.94
end
"""
    approximate_surface_area(a::Real, R::Real ,κ::Real, δ::Real)

Approximation of the surface area of a miller geometry flux surface
"""
function approximate_surface_area(a::Real, R::Real, κ::Real, δ::Real)
    return 2pi^2 * R * a^2 * κ * (1.0 - 0.151 * δ * a / R) / (2pi * R) # m²
end

"""
    silo(r_start, r_end, height_start, height_end; n_points::Int=100, resolution::Float64=1.0)
"""
function silo(r_start::Real, r_end::Real, height_start::Real, height_end::Real; n_points::Int=100, resolution::Float64=1.0)
    n_points = Int(round(n_points * resolution))
    height_start = abs(height_start)
    height_end = abs(height_end)
    height_end = min(max(height_end, height_start * 0.0), height_start * 0.9)
    x, y = ellipse(r_end - r_start, height_start - height_end, 0.0, pi / 2, r_start, height_end; n_points)
    return vcat(r_start, r_start, r_end, x), vcat(height_start, 0.0, 0.0, y) .- (height_start / 2.0)
end

"""
    line_through_point(m, x0, y0, x)

return `y` values at `x` of line with gradient `m` going through point `(x0,y0)`
"""
function line_through_point(m::T, x0::T, y0::T, x::Vector{T}) where {T<:Real}
    return @. m * x + y0 - m * x0
end

mutable struct MXHboundary
    mxh::IMAS.MXH
    RX::Vector{<:Real}
    ZX::Vector{<:Real}
    r_boundary::Vector{<:Real}
    z_boundary::Vector{<:Real}
end

@recipe function plot_mxhb(mxhb::MXHboundary)
    @series begin
        aspect_ratio --> :equal
        label --> ""
        lw --> 2
        color --> :black
        mxhb.r_boundary, mxhb.z_boundary
    end
    @series begin
        primary := false
        aspect_ratio --> :equal
        label --> ""
        mxhb.mxh
    end
    @series begin
        seriestype := :scatter
        primary := false
        aspect_ratio --> :equal
        label --> ""
        marker := :cross
        mxhb.RX, mxhb.ZX
    end
end

"""
    add_xpoint(
        mr::AbstractVector{T},
        mz::AbstractVector{T},
        R0::Union{Nothing,T}=nothing,
        Z0::Union{Nothing,T}=nothing;
        upper::Bool,
        α_multiplier::Float64
    ) where {T<:Real}

Add a X-point to a boundary that does not have one.

The X-point is added at the point where there is the largest curvature in the boundary,
and is built in such a way to form 90° angle (which is what it should be for zero current at the LCFS).
The X-point is placed on a line that goes from (R0,Z0) through the point of maximum curvature.

Control of the X-point location can be achieved by modifying R0, Z0.
"""
function add_xpoint(
    mr::AbstractVector{T},
    mz::AbstractVector{T},
    R0::T,
    Z0::T;
    upper::Bool,
    α_multiplier::Float64
) where {T<:Real}

    function cost(pr::Vector{Float64}, pz::Vector{Float64}, i::Integer, R0::T, Z0::T, α::Float64)
        return abs(add_xpoint(pr, pz, i, R0, Z0, α).θX - π / 2)
    end

    if upper
        i = argmax(abs.(IMAS.curvature(mr, mz)) .* (mz .> Z0))
    else
        i = argmax(abs.(IMAS.curvature(mr, mz)) .* (mz .< Z0))
    end

    res = Optim.optimize(α -> cost(mr, mz, i, R0, Z0, α), 1.0, 1.5, Optim.GoldenSection())
    return add_xpoint(mr, mz, i, R0, Z0, 1.0 + (res.minimizer[1] - 1.0) * α_multiplier)
end

function add_xpoint(pr::Vector{Float64}, pz::Vector{Float64}, i::Integer, R0::T, Z0::T, α::T) where {T<:Real}
    @assert IMAS.is_closed_polygon(pr, pz)
    RX = pr[i] .* α .+ R0 .* (1.0 .- α)
    ZX = pz[i] .* α .+ Z0 .* (1.0 .- α)

    pr = @views pr[1:end-1]
    pz = @views pz[1:end-1]
    n = length(pr)

    min_angle = Inf
    i_min = 0
    max_angle = -Inf
    i_max = 0
    for i in 1:n
        @inbounds angle = atan(pz[i] - ZX, pr[i] - RX)
        if angle < min_angle
            min_angle = angle
            i_min = i
        end
        if angle > max_angle
            max_angle = angle
            i_max = i
        end
    end

    index = [mod(idx - 1, n) + 1 for idx in i_max:i_max+n-(i_max-i_min)]
    R = @views [RX; pr[index]; RX]
    Z = @views [ZX; pz[index]; ZX]

    θX = max_angle - min_angle

    return (RX=RX, ZX=ZX, R=R, Z=Z, θX=θX)
end

"""
    MXHboundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, n_points::Integer=0)

Return MXHboundary structure of boundary with x-points based on input MXH boundary parametrization
"""
function MXHboundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, n_points::Integer=0)
    mr, mz = mxh()
    mxhb = MXHboundary(deepcopy(mxh), Float64[], Float64[], mr, mz)
    return MXHboundary!(mxhb; upper_x_point, lower_x_point, n_points)
end

function MXHboundary!(mxhb::MXHboundary; upper_x_point::Bool, lower_x_point::Bool, n_points::Integer=0, RU::Real=0.0, RL::Real=0.0)
    mr, mz = mxhb.mxh()

    if ~upper_x_point && ~lower_x_point && n_points === 0
        empty!(mxhb.RX)
        empty!(mxhb.ZX)
        mxhb.r_boundary = mr
        mxhb.z_boundary = mz
        return mxhb
    end

    R0 = mxhb.mxh.R0
    Z0 = mxhb.mxh.Z0
    if RU == 0.0
        RU = R0
    end
    if RL == 0.0
        RL = R0
    end

    R, Z = mr, mz
    RXU, ZXU, R1, Z1, _ = add_xpoint(R, Z, RU, Z0; upper=true, α_multiplier=(upper_x_point ? 1.0 : 2.0))
    if upper_x_point
        R = R1
        Z = Z1
    end
    RXL, ZXL, R2, Z2 = add_xpoint(R, Z, RL, Z0; upper=false, α_multiplier=(lower_x_point ? 1.0 : 2.0))
    if lower_x_point
        R = R2
        Z = Z2
    end

    RX = Float64[]
    ZX = Float64[]
    if upper_x_point
        push!(RX, RXU)
        push!(ZX, ZXU)
    end
    if lower_x_point
        push!(RX, RXL)
        push!(ZX, ZXL)
    end

    # resample boundary after convex_hull in such a way to preserve x-points and add proper curvature in the x-point region
    if upper_x_point + lower_x_point == 1
        if upper_x_point
            R, Z = IMAS.reorder_flux_surface!(R, Z, argmax(Z))
        else
            R, Z = IMAS.reorder_flux_surface!(R, Z, argmin(Z))
        end
        RR = @views [-(reverse(R[2:end-1]) .- R[1]) .+ R[1]; R[2:end-1]; -(reverse(R[2:end-1]) .- R[1]) .+ R[1]]
        ZZ = @views [-(reverse(Z[2:end-1]) .- Z[1]) .+ Z[1]; Z[2:end-1]; -(reverse(Z[2:end-1]) .- Z[1]) .+ Z[1]]
        RR, ZZ = IMAS.resample_plasma_boundary(RR, ZZ; n_points=length(R) * 2, method=:cubic)
        if upper_x_point
            I = ZZ .< Z[1]
        else
            I = ZZ .> Z[1]
        end
        R = [R[1]; RR[I]; R[1]]
        Z = [Z[1]; ZZ[I]; Z[1]]
    else
        R, Z = IMAS.reorder_flux_surface!(R, Z, argmax(Z))
        izmin = argmin(Z)
        R1 = @views R[1:izmin]
        Z1 = @views Z[1:izmin]
        R2 = @views R[izmin:end]
        Z2 = @views Z[izmin:end]
        ΔZ = (Z1[end] - Z[1])
        ZZ1 = @views [Z1[2:end-1] .- ΔZ; Z1[2:end-1]; Z1[2:end-1] .+ ΔZ]
        RR1 = @views [-(reverse(R1[2:end-1]) .- R1[1]) .+ R1[1]; R1[2:end-1]; -(reverse(R1[2:end-1]) .- R1[1]) .+ R1[1]]
        ZZ2 = @views [Z2[2:end-1] .+ ΔZ; Z2[2:end-1]; Z2[2:end-1] .- ΔZ]
        RR2 = @views [-(reverse(R2[2:end-1]) .- R2[1]) .+ R2[1]; R2[2:end-1]; -(reverse(R2[2:end-1]) .- R2[1]) .+ R2[1]]
        RR1, ZZ1 = IMAS.resample_2d_path(RR1, ZZ1; n_points=length(R) * 3, method=:cubic)
        RR2, ZZ2 = IMAS.resample_2d_path(RR2, ZZ2; n_points=length(R) * 3, method=:cubic)
        crossings = IMAS.intersection(RR1, ZZ1, RR2, ZZ2).crossings
        if crossings[1][2] > crossings[2][2]
            ((RXU, ZXU), (RXL, ZXL)) = crossings
        else
            ((RXU, ZXU), (RXL, ZXL)) = crossings
        end
        I1 = (ZZ1 .< ZXU) .&& (ZZ1 .> ZXL)
        I2 = (ZZ2 .< ZXU) .&& (ZZ2 .> ZXL)
        R = [RXU; RR1[I1]; RXL; RR2[I2]; RXU]
        Z = [ZXU; ZZ1[I1]; ZXL; ZZ2[I2]; ZXU]
    end

    IMAS.reorder_flux_surface!(R, Z, R0, Z0)

    mxhb.RX = [RXU, RXL]
    mxhb.ZX = [ZXU, ZXL]
    mxhb.r_boundary = R
    mxhb.z_boundary = Z

    return mxhb
end

function fitMXHboundary(mxh::IMAS.MXH, nx::Int; n_points::Int=0)
    return fitMXHboundary(mxh; upper_x_point=nx ∈ (1, 2), lower_x_point=nx ∈ (-1, 2), n_points)
end

"""
    fitMXHboundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, n_points::Integer=0)

Find boundary such that the output MXH parametrization (with x-points) matches the input MXH parametrization (without x-points)
"""
function fitMXHboundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, n_points::Int=0, target_area::Float64=0.0, target_volume::Float64=0.0, debug::Bool=false)
    if ~upper_x_point && ~lower_x_point
        return MXHboundary(mxh; upper_x_point, lower_x_point, n_points)
    end

    symmetric = IMAS.is_updown_symmetric(mxh) && !xor(upper_x_point, lower_x_point)

    pr, pz = mxh()
    i = argmax(pz)
    RXU = pr[i]
    ZXU = pz[i]
    i = argmin(pz)
    RXL = pr[i]
    ZXL = pz[i]
    pr, pz = IMAS.resample_plasma_boundary(pr, pz; n_points=100)

    M = 12
    N = min(M, length(mxh.s))

    # change mxhb0 parameters so that the boundary  with x-points fit has the desired elongation, triangularity, squareness when fit with mxh0 
    mxhb0 = MXHboundary(mxh; upper_x_point, lower_x_point, n_points)
    mxh0 = IMAS.MXH(mxhb0.r_boundary, mxhb0.z_boundary, M)

    function mxhb_from_params!(mxhb0::MXHboundary, params::AbstractVector{<:Real}; upper_x_point::Bool, lower_x_point::Bool, n_points::Int)
        L = 1
        mxhb0.mxh.Z0 = params[L]
        L += 1
        mxhb0.mxh.κ = params[L]
        if symmetric
            mxhb0.mxh.c0 = 0.0
            mxhb0.mxh.s = params[L+1:end]
            mxhb0.mxh.c = mxhb0.mxh.s .* 0.0
        else
            L += 1
            mxhb0.mxh.c0 = params[L]
            mxhb0.mxh.s = params[L+1:L+Integer((end - L) / 2)]
            mxhb0.mxh.c = params[L+Integer((end - L) / 2)+1:end]
        end
        return MXHboundary!(mxhb0; upper_x_point, lower_x_point, n_points)
    end

    function cost(params::AbstractVector{<:Real}; mxhb0, mxh0, upper_x_point, lower_x_point, n_points, target_area::Float64, target_volume::Float64)
        # mxhb0: contains the boundary with the x-point
        # mxhb0.mxh: is the MXH parametrization that leads to mxhb0 once x-points are set
        # mxh0: is the MHX fit to the mxhb0 boundary

        mxhb_from_params!(mxhb0, params; upper_x_point, lower_x_point, n_points)
        IMAS.MXH!(mxh0, mxhb0.r_boundary, mxhb0.z_boundary)

        pr0, pz0 = IMAS.resample_plasma_boundary(mxhb0.r_boundary, mxhb0.z_boundary; n_points=length(pr))

        # X-points and non X-points halves
        c = 0.0
        if upper_x_point
            i = argmax(mxhb0.ZX)
            c += ((mxhb0.RX[i] - RXU)^2 + (mxhb0.ZX[i] - ZXU)^2) / mxhb0.mxh.R0^2
        else
            i = pz .< mxhb0.mxh.Z0
            c += sum((pr0[i] .- pr[i]) .^ 2 .+ (pz0[i] .- pz[i]) .^ 2) / mxhb0.mxh.R0^2 / sum(i)
        end
        if lower_x_point
            i = argmin(mxhb0.ZX)
            c += ((mxhb0.RX[i] - RXL)^2 + (mxhb0.ZX[i] - ZXL)^2) / mxhb0.mxh.R0^2
        else
            i = pz .> mxhb0.mxh.Z0
            c += sum((pr0[i] .- pr[i]) .^ 2 .+ (pz0[i] .- pz[i]) .^ 2) / mxhb0.mxh.R0^2 / sum(i)
        end

        # other shape parameters
        c += (mxh0.R0 - mxh.R0)^2
        c += (mxh0.Z0 - mxh.Z0)^2
        c += (mxh0.κ - mxh.κ)^2
        c += (mxh0.ϵ - mxh.ϵ)^2
        c += (mxh0.c0 - mxh.c0)^2
        c += sum((mxh0.s[1:N] .- mxh.s[1:N]) .^ 2)
        c += sum((mxh0.c[1:N] .- mxh.c[1:N]) .^ 2)

        # make MXH match boundary of MHXboundary
        x_point_area = IMAS.area(mxhb0.r_boundary, mxhb0.z_boundary)
        marea0 = IMAS.area(mxh0()...)
        c += (abs(x_point_area - marea0) / x_point_area)^2 * 1E3

        # To avoid MXH solutions with kinks force area and convex_hull area to match
        mr, mz = mxhb0.mxh()
        hull = convex_hull(mr, mz; closed_polygon=true)
        mrch = [r for (r, z) in hull]
        mzch = [z for (r, z) in hull]
        marea = IMAS.area(mr, mz)
        mareach = IMAS.area(mrch, mzch)
        c += (abs(mareach - marea) / mareach)^2 * 1E3

        # Matching a target area and or volume
        if target_area > 0
            c += ((x_point_area - target_area)^2 * 1E2)
        end
        if target_volume > 0
            x_point_volume = IMAS.revolution_volume(mxhb0.r_boundary, mxhb0.z_boundary)
            c += ((x_point_volume - target_volume)^2 * 1E1)
        end

        return sqrt(c)
    end

    if symmetric
        params = vcat(mxh.Z0, mxh.κ, mxh.s)
    else
        params = vcat(mxh.Z0, mxh.κ, mxh.c0, mxh.s, mxh.c)
    end
    res = Optim.optimize(
        x -> cost(x; mxhb0, mxh0, upper_x_point, lower_x_point, n_points, target_area, target_volume),
        params,
        Optim.NelderMead(),
        Optim.Options(; iterations=1000, g_tol=1E-5)
    )
    mxhb_from_params!(mxhb0, res.minimizer; upper_x_point, lower_x_point, n_points)
    IMAS.MXH!(mxh0, mxhb0.r_boundary, mxhb0.z_boundary)

    if debug
        println(res)

        println("mxh")
        println(mxh)

        println("mxh0")
        println(mxh0)

        println("mxhb0.mxh")
        println(mxhb0.mxh)

        display(plot(mxh0))
    end

    return mxhb0
end

"""
    free_boundary_private_flux_constraint(R::AbstractVector{T}, Z::AbstractVector{T}; upper_x_point::Bool, lower_x_point::Bool, fraction::Float64=0.25, n_points::Int=10) where {T<:Real}

Given a closed boundary R,Z returns Rp,Zp constraint points on plausible private flux regions
"""
function free_boundary_private_flux_constraint(
    R::AbstractVector{T},
    Z::AbstractVector{T};
    upper_x_point::Bool,
    lower_x_point::Bool,
    fraction::Float64=0.25,
    n_points::Int=10) where {T<:Real}

    Rp = T[]
    Zp = T[]

    if fraction <= 0.0
        return Rp, Zp
    end

    for x_point in [-1, 1]

        if x_point == -1 && !lower_x_point
            continue
        elseif x_point == 1 && !upper_x_point
            continue
        end

        pr = deepcopy(R)
        pz = deepcopy(Z)

        # up-down and left-right mirror of lcfs with respect to the X-point
        if x_point == 1
            pr, pz = IMAS.reorder_flux_surface!(pr, pz, argmin(pz))
            Rx = pr[argmax(pz)]
            Zx = pz[argmax(pz)]
        else
            pr, pz = IMAS.reorder_flux_surface!(pr, pz, argmax(pz))
            Rx = pr[argmin(pz)]
            Zx = pz[argmin(pz)]
        end
        pr = -(pr .- Rx) .+ Rx
        pz = -(pz .- Zx) .+ Zx

        # select points withing a given radius
        d = sqrt.((pr .- Rx) .^ 2 .+ (pz .- Zx) .^ 2)
        index = d .< abs(maximum(pz) - minimum(pz)) * fraction
        pr = pr[index]
        pz = pz[index]

        # resample to give requested number of points
        pr, pz = IMAS.resample_plasma_boundary(pr, pz; n_points)

        append!(Rp, pr)
        append!(Zp, pz)
    end

    return Rp, Zp
end

function boundary_shape(R0::Real; p=nothing)
    return boundary_shape(IMAS.MXH(R0, 3); p=p)
end

"""
    private_flux_regions_from_lcfs(mr::AbstractArray{T}, mz::AbstractArray{T}, len::T; n_points::Integer=2)

Traces private flux regions starting from LCFS boundary
"""
function private_flux_regions_from_lcfs(mr::AbstractArray{T}, mz::AbstractArray{T}, len::T; n_points::Integer=2) where {T<:Real}
    mr = deepcopy(mr)
    mz = deepcopy(mz)
    IMAS.reorder_flux_surface!(mr, mz, sum(mr) / length(mr), sum(mz) / length(mz))

    iu = argmax(mz)
    il = argmin(mz)

    r_hfs = reverse(mr[iu:il])
    z_hfs = reverse(mz[iu:il])

    r_lfs = vcat(mr[il:end-1], mr[1:iu])
    z_lfs = vcat(mz[il:end-1], mz[1:iu])

    zz_u_hfs = range(mz[iu] + len, mz[iu], n_points)
    zz_u_lfs = range(mz[iu], mz[iu] + len, n_points)
    zz_l_hfs = range(mz[il] - len, mz[il], n_points)
    zz_l_lfs = range(mz[il], mz[il] - len, n_points)
    rr_u_hfs = IMAS.interp1d(z_lfs, r_lfs).(zz_u_hfs)
    rr_u_lfs = IMAS.interp1d(z_hfs, r_hfs).(zz_u_lfs)
    rr_l_hfs = IMAS.interp1d(z_lfs, r_lfs).(zz_l_hfs)
    rr_l_lfs = IMAS.interp1d(z_hfs, r_hfs).(zz_l_lfs)

    rr_u = vcat(rr_u_hfs[1:end-1], rr_u_lfs[2:end])
    zz_u = vcat(zz_u_hfs[1:end-1], zz_u_lfs[2:end])
    rr_l = vcat(rr_l_hfs[1:end-1], rr_l_lfs[2:end])
    zz_l = vcat(zz_l_hfs[1:end-1], zz_l_lfs[2:end])

    return rr_u, zz_u, rr_l, zz_l
end

"""
    check_ped_finder(profile::AbstractVector{<:Real},psi_norm::AbstractVector{<:Real})

Plots the profile and fitted profile for pedestal finding (region outside pedestal not important for fit)
"""
function check_ped_finder(profile::AbstractVector{<:Real}, psi_norm::AbstractVector{<:Real})
    ped_height, ped_width = IMAS.pedestal_finder(profile, psi_norm)
    plot(psi_norm, profile; label="original profile")
    return plot!(
        psi_norm,
        IMAS.Hmode_profiles(profile[end], ped_height, profile[1], length(profile), 2.0, 2.0, ped_width);
        label="fitted profile (pedestal region is important only)"
    )
end

"""
    ip_from_q_star(a::Real, B0::Real, q_star::Real, R0::Real, κ::Real)

Calculate ip from the cylindrical equivalent safety factor
"""
function ip_from_q_cyl(a::Real, B0::Real, q_cyl::Real, R0::Real, κ::Real)
    return 1e6 * abs((5 * a^2 * B0) / (R0 * (q_cyl)) * (1 + κ^2) / 2)
end

"""
    B0_from_Bmax(Bmax::Real, distance_tf_to_plasma::Real, a::Real, R0::Real)

Calculate the on axis magnetic field from the maximum magnetic field at the tf coils
"""
function B0_from_b_max(b_max::Real, distance_tf_to_plasma::Real, a::Real, R0::Real)
    return b_max * (1 - (a + distance_tf_to_plasma) / R0)
end