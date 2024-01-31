import SpecialFunctions
import QuadGK
import Interpolations
import PolygonOps

#= =============== =#
#  Shape functions  #
#= =============== =#

function layer_shape_message(shape_function_index)
    return "layer.shape=$(shape_function_index) is invalid. Valid options are:
            -2: Offset & convex-hull
            -1: Offset
             1: Priceton D exact  (shape_parameters = [])
             2: Priceton D approx (shape_parameters = [])
             3: Priceton D scaled (shape_parameters = [height])
             4: rectangle         (shape_parameters = [height])
             5: double_ellipse    (shape_parameters = [centerpost_height, height])
             6: tripple-arc       (shape_parameters = [height, small_radius, mid_radius, small_coverage, mid_coverage])
             7: miller            (shape_parameters = [elongation, triangularity])
             8: square_miller     (shape_parameters = [elongation, triangularity, squareness])
             9: spline            (shape_parameters = [hfact, rz...)
            10: silo              (shape_parameters = [h_start, h_end)
           10x: shape + z_offset  (shape_parameters = [..., z_offset])
          100x: negative shape    (shape_parameters = [...])"
end

function initialize_shape_parameters(shape_function_index, r_obstruction, z_obstruction, r_start, r_end, target_clearance)
    height = maximum(z_obstruction) - minimum(z_obstruction) + target_clearance * 2.0
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
            centerpost_height = (maximum(z_obstruction) - minimum(z_obstruction))
            shape_parameters = [centerpost_height, height]
        elseif shape_index_mod == Int(_rectangle_)
            shape_parameters = [height]
        elseif shape_index_mod == Int(_triple_arc_)
            shape_parameters = [log10(height), log10(1E-3), log10(1E-3), log10(45), log10(45)]
        elseif shape_index_mod ∈ (Int(_miller_), Int(_square_miller_))
            _, imaxr = findmax(r_obstruction)
            _, iminr = findmin(r_obstruction)
            _, imaxz = findmax(z_obstruction)
            _, iminz = findmin(z_obstruction)
            r_at_max_z, max_z = r_obstruction[imaxz], z_obstruction[imaxz]
            r_at_min_z, min_z = r_obstruction[iminz], z_obstruction[iminz]
            z_at_max_r, max_r = z_obstruction[imaxr], r_obstruction[imaxr]
            z_at_min_r, min_r = z_obstruction[iminr], r_obstruction[iminr]
            a = 0.5 * (max_r - min_r)
            b = 0.5 * (max_z - min_z)
            R = 0.5 * (max_r + min_r)
            elongation = b / a
            triup = (R - r_at_max_z) / a
            tridown = (R - r_at_min_z) / a
            if shape_index_mod == Int(_miller_)
                shape_parameters = [elongation, (triup + tridown) / 2.0]
            else
                shape_parameters = [elongation, (triup + tridown) / 2.0, 0.0]
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
        elseif shape_index_mod == Int(_double_ellipse_)
            func = circle_ellipse
        elseif shape_index_mod == Int(_rectangle_)
            func = rectangle_shape
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
    function resampled_zfunc(args...; resample=true)
        if resample
            return IMAS.resample_2d_path(dfunc(args...)...; method=:linear)
        else
            return dfunc(args...)
        end
    end

    return resampled_zfunc
end

"""
    optimize_shape(r_obstruction, z_obstruction, target_clearance, func, r_start, r_end, shape_parameters; verbose=false)

Find shape parameters that generate smallest shape and target clearance from an obstruction
"""
function optimize_shape(
    r_obstruction::Vector{Float64},
    z_obstruction::Vector{Float64},
    target_clearance::Float64,
    func::Function,
    r_start::Float64,
    r_end::Float64,
    shape_parameters::Vector{Float64};
    use_curvature::Bool=true,
    verbose::Bool=false
)

    rz_obstruction = collect(zip(r_obstruction, z_obstruction))
    initial_guess = deepcopy(shape_parameters)

    if length(shape_parameters) in (0, 1)
        func(r_start, r_end, shape_parameters...)

    else
        function cost_shape(
            r_obstruction::Vector{Float64},
            z_obstruction::Vector{Float64},
            rz_obstruction::Vector{Tuple{Float64,Float64}},
            target_clearance::Float64,
            target_area::Float64,
            func::Function,
            r_start::Float64,
            r_end::Float64,
            shape_parameters::Vector{Float64};
            use_curvature::Bool,
            verbose::Bool=false
        )
            R, Z = func(r_start, r_end, shape_parameters...)

            cost_area = abs(IMAS.area(R, Z) - target_area) / target_area

            length(R) == length(Z) || error("R and Z have different length")
            # disregard near r_start and r_end where optimizer has no control and shape is allowed to go over obstruction
            index = (.!)(isapprox.(R, r_start) .|| isapprox.(R, r_end))
            Rv = view(R, index)
            Zv = view(Z, index)

            # no polygon crossings  O(N)
            cost_inside = 0.0
            for (r, z) in zip(Rv, Zv)
                inpoly = PolygonOps.inpolygon((r, z), rz_obstruction)
                cost_inside += inpoly
            end

            # target clearance  O(1)
            minimum_distance, mean_above_distance_error = IMAS.min_mean_distance_error_two_shapes(Rv, Zv, r_obstruction, z_obstruction, target_clearance; above_target=true)
            if minimum_distance < target_clearance
                cost_min_clearance = (minimum_distance - target_clearance) / target_clearance
            else
                cost_min_clearance = 0.0
            end
            cost_mean_above_distance = mean_above_distance_error / target_clearance

            # curvature
            if use_curvature
                curvature = abs.(IMAS.curvature(Rv, Zv))
                cost_max_curvature = exp(maximum(curvature)^2.5)
            else
                cost_max_curvature = 0.0
            end

            # favor up/down symmetric solutions
            cost_up_down_symmetry = abs(maximum(Z) + minimum(Z)) / (maximum(Z) - minimum(Z))

            if verbose
                @show minimum_distance
                @show mean_above_distance_error
                @show target_clearance
                @show cost_min_clearance^2
                @show cost_mean_above_distance^2
                @show cost_inside^2
                @show cost_up_down_symmetry^2
                @show cost_max_curvature^2
            end

            # return cost
            return cost_min_clearance^2 + cost_mean_above_distance^2 + cost_inside^2 + 0.1 * cost_up_down_symmetry^2 + cost_area + cost_max_curvature
        end

        r_obstruction_buffered, z_obstruction_buffered = buffer(r_obstruction, z_obstruction, target_clearance)
        target_area = IMAS.area(r_obstruction_buffered, z_obstruction_buffered)

        initial_guess = copy(shape_parameters)
        res = Optim.optimize(
            shape_parameters -> cost_shape(r_obstruction, z_obstruction, rz_obstruction, target_clearance, target_area, func, r_start, r_end, shape_parameters; use_curvature),
            initial_guess, length(shape_parameters) == 1 ? Optim.BFGS() : Optim.NelderMead(), Optim.Options())
        if verbose
            println(res)
        end
        shape_parameters = Optim.minimizer(res)
    end

    # check no polygon crossings
    R, Z = func(r_start, r_end, shape_parameters...)
    cost_inside = 0
    for (r, z) in zip(R, Z)
        inpoly = PolygonOps.inpolygon((r, z), rz_obstruction)
        cost_inside += inpoly
    end

    if cost_inside > 0
        @warn "optimize_shape function could not avoid polygon crossings! Perhaps try changing shape?"
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
    return double_ellipse(r_start, r_end, r_center, centerpost_height, 0.0, height; n_points, resolution)
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
    centerpost_height = abs(centerpost_height)
    height = abs(height)
    centerpost_height = mirror_bound(centerpost_height, 0.0, height)
    r_center = r_start + (height - centerpost_height) / 2.0
    return double_ellipse(r_start, r_end, r_center, centerpost_height, height; n_points, resolution)
end

"""
    rectangle_shape(r_start::T, r_end::T, z_low::T, z_high::T; n_points::Int=5) where {T<:Real}

Asymmetric rectangular shape
"""
function rectangle_shape(r_start::T, r_end::T, z_low::T, z_high::T; n_points::Int=5, resolution::Float64=1.0) where {T<:Real}
    n_points = Int(round(n_points * resolution))
    n_points = max(5, n_points)
    if n_points == 5
        R = [r_start, r_start, r_end, r_end, r_start]
        Z = [z_low, z_high, z_high, z_low, z_low]
    else
        R = vcat(
            r_start,
            range(r_start, r_start, n_points)[2:end],
            range(r_end, r_start, n_points)[2:end],
            range(r_end, r_end, n_points)[2:end],
            range(r_start, r_end, n_points))
        Z = vcat(
            z_low,
            range(z_high, z_low, n_points)[2:end],
            range(z_high, z_high, n_points)[2:end],
            range(z_low, z_high, n_points)[2:end],
            range(z_low, z_low, n_points))
    end
    return R, Z
end

"""
    rectangle_shape(r_start::T, r_end::T, height::T; n_points::Integer=5) where {T<:Real}

Symmetric rectangular shape
"""
function rectangle_shape(r_start::T, r_end::T, height::T; n_points::Integer=5, resolution::Float64=1.0) where {T<:Real}
    Δ = height / 2.0
    return rectangle_shape(r_start, r_end, -Δ, Δ; n_points, resolution)
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

height, small_radius, mid_radius, small_coverage, mid_coverage are 10^exponent (to ensure positiveness)
"""
function triple_arc(
    r_start::T,
    r_end::T,
    height::T,
    small_radius::T,
    mid_radius::T,
    small_coverage::T,
    mid_coverage::T;
    min_small_radius_fraction::T=0.2,
    min_mid_radius_fraction::T=min_small_radius_fraction * 2.0,
    n_points::Integer=400,
    resolution::Float64=1.0) where {T<:Real}

    n_points = Int(round(n_points / 4.0 * resolution))

    height = 10^height / 2.0
    small_radius = 10^small_radius + height * min_small_radius_fraction
    mid_radius = 10^mid_radius + height * min_mid_radius_fraction
    small_coverage = 10^small_coverage * pi / 180
    mid_coverage = 10^mid_coverage * pi / 180

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
    n_points = Int(floor(n_points * resolution) / 2) * 2 + 1

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

    n_points = Int(floor(n_points * resolution) / 2) * 2 + 1

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
    x_b = T[v[1] for v in GeoInterface.coordinates(poly_b)[1]]
    y_b = T[v[2] for v in GeoInterface.coordinates(poly_b)[1]]
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
    x_b = T[v[1] for v in GeoInterface.coordinates(poly_b)[1]]
    y_b = T[v[2] for v in GeoInterface.coordinates(poly_b)[1]]
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
    add_xpoint(mr::Vector{T}, mz::Vector{T}, R0::Union{Nothing,T}, Z0::T; upper::Bool) where {T<:Real}

Add a X-point to a boundary that does not have one.

The X-point is added at the point where there is the largest curvature in the boundary,
and is built in such a way to form 90° angle (which is what it should be for zero current at the LCFS).
The X-point is placed on a line that goes from (R0,Z0) through the point of maximum curvature.

Control of the X-point location can be achieved by modifying R0, Z0.
"""
function add_xpoint(mr::AbstractVector{T}, mz::AbstractVector{T}, R0::Union{Nothing,T}=nothing, Z0::Union{Nothing,T}=nothing; upper::Bool) where {T<:Real}

    # NOTE: when calling curvature we multiply Z by 2 to form x-points preferentially in Z direction (even for circular plasma, for example)
    function cost(points0::Vector, points::Vector, i::Integer, R0::T, Z0::T, α::Float64)
        points .= points0
        if upper
            RX, ZX, R, Z = add_xpoint(points, i, R0, Z0, α)
            return (1.0 - maximum(abs, IMAS.curvature(R, 2 * Z)[(Z.>(Z0+ZX)/2.0)]))^2.0
        else
            RX, ZX, R, Z = add_xpoint(points, i, R0, Z0, α)
            return (1.0 - maximum(abs, IMAS.curvature(R, 2 * Z)[(Z.<(Z0+ZX)/2.0)]))^2.0
        end
    end

    if Z0 === nothing
        Z0 = sum(mz) / length(mz)
    end

    if upper
        k = argmax(abs.(IMAS.curvature(mr, 2 * mz)) .* (mz .> Z0))
        index = mz .> Z0
        mri = mr[index]
        mzi = mz[index]
        i = findfirst(i -> (mri[i] == mr[k] && mzi[i] == mz[k]), 1:length(mri))
    else
        k = argmax(abs.(IMAS.curvature(mr, 2 * mz)) .* (mz .< Z0))
        index = mz .< Z0
        mri = mr[index]
        mzi = mz[index]
        i = findfirst(i -> (mri[i] == mr[k] && mzi[i] == mz[k]), 1:length(mri))
    end

    if R0 === nothing
        R0 = mri[i]
    end

    points = collect(zip([mri; mri[1]], [mzi; mzi[1]]))
    points0 = deepcopy(points)

    res = Optim.optimize(α -> cost(points0, points, i, R0, Z0, α), 1.05, 1.5, Optim.GoldenSection())

    points = collect(zip([mr; mr[1]], [mz; mz[1]]))
    RX, ZX, R, Z = add_xpoint(points, k, R0, Z0, res.minimizer[1])

    return RX, ZX, R, Z
end

function add_xpoint(points::Vector, i::Integer, R0::T, Z0::T, α::T) where {T<:Real}
    RX = points[i][1] .* α .+ R0 .* (1.0 .- α)
    ZX = points[i][2] .* α .+ Z0 .* (1.0 .- α)
    points[end] = (RX, ZX)
    RZ = convex_hull!(points; closed_polygon=true)
    R = T[r for (r, z) in RZ]
    Z = T[z for (r, z) in RZ]
    return RX, ZX, R, Z
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

    RX = Float64[]
    ZX = Float64[]
    if upper_x_point
        RXU, ZXU, _ = add_xpoint(mr, mz, RU, Z0; upper=true)
        push!(RX, RXU)
        push!(ZX, ZXU)
    end
    if lower_x_point
        RXL, ZXL, _ = add_xpoint(mr, mz, RL, Z0; upper=false)
        push!(RX, RXL)
        push!(ZX, ZXL)
    end

    RZ = convex_hull(vcat(mr, RX), vcat(mz, ZX); closed_polygon=true)
    R = [r for (r, z) in RZ]
    Z = [z for (r, z) in RZ]

    # resample boundary after convex_hull in such a way to preserve x-points and add proper curvature in the x-point region
    if upper_x_point + lower_x_point == 1
        if upper_x_point
            R, Z = IMAS.reorder_flux_surface!(R, Z, argmax(Z))
        else
            R, Z = IMAS.reorder_flux_surface!(R, Z, argmin(Z))
        end
        RR = [-(reverse(R[2:end-1]) .- R[1]) .+ R[1]; R[2:end-1]; -(reverse(R[2:end-1]) .- R[1]) .+ R[1]]
        ZZ = [-(reverse(Z[2:end-1]) .- Z[1]) .+ Z[1]; Z[2:end-1]; -(reverse(Z[2:end-1]) .- Z[1]) .+ Z[1]]
        RR, ZZ = IMAS.resample_plasma_boundary(RR, ZZ; n_points=length(R) * 2)
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
        R1, R2 = R[1:izmin], R[izmin:end]
        Z1, Z2 = Z[1:izmin], Z[izmin:end]
        ΔZ = (Z1[end] - Z[1])
        ZZ1 = [Z1[2:end-1] .- ΔZ; Z1[2:end-1]; Z1[2:end-1] .+ ΔZ]
        RR1 = [-(reverse(R1[2:end-1]) .- R1[1]) .+ R1[1]; R1[2:end-1]; -(reverse(R1[2:end-1]) .- R1[1]) .+ R1[1]]
        ZZ2 = [Z2[2:end-1] .+ ΔZ; Z2[2:end-1]; Z2[2:end-1] .- ΔZ]
        RR2 = [-(reverse(R2[2:end-1]) .- R2[1]) .+ R2[1]; R2[2:end-1]; -(reverse(R2[2:end-1]) .- R2[1]) .+ R2[1]]
        RR1, ZZ1 = IMAS.resample_plasma_boundary(RR1, ZZ1; n_points=length(R) * 2)
        RR2, ZZ2 = IMAS.resample_plasma_boundary(RR2, ZZ2; n_points=length(R) * 2)
        I1 = (ZZ1 .< Z1[1]) .&& (ZZ1 .> Z1[end])
        I2 = (ZZ2 .< Z1[1]) .&& (ZZ2 .> Z1[end])
        R = [R1[1]; RR1[I1]; R1[end]; RR2[I2]; R1[1]]
        Z = [Z1[1]; ZZ1[I1]; Z1[end]; ZZ2[I2]; Z1[1]]
    end

    IMAS.reorder_flux_surface!(R, Z, R0, Z0)

    mxhb.RX = RX
    mxhb.ZX = ZX
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

    # change mxh parameters so that the boundary with x-points has the desired elongation, triangularity, squareness
    mxhb0 = MXHboundary(mxh; upper_x_point, lower_x_point, n_points)
    mxh0 = IMAS.MXH(mxhb0.r_boundary, mxhb0.z_boundary, M)

    function mxhb_from_params!(mxhb0::MXHboundary, params::AbstractVector{<:Real}; upper_x_point::Bool, lower_x_point::Bool, n_points::Int)
        L = 1
        mxhb0.mxh.Z0 = params[L]
        L += 1
        mxhb0.mxh.κ = params[L]
        L += 1
        mxhb0.mxh.c0 = params[L]
        mxhb0.mxh.s = params[L+1:L+Integer((end - L) / 2)]
        mxhb0.mxh.c = params[L+Integer((end - L) / 2)+1:end]
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
    res = Optim.optimize(
        x -> cost(x; mxhb0, mxh0, upper_x_point, lower_x_point, n_points, target_area, target_volume),
        vcat(mxh.Z0, mxh.κ, mxh.c0, mxh.s, mxh.c),
        Optim.NelderMead(),
        Optim.Options(; iterations=1000)
    )
    mxhb_from_params!(mxhb0, res.minimizer; upper_x_point, lower_x_point, n_points)
    IMAS.MXH!(mxh0, mxhb0.r_boundary, mxhb0.z_boundary)
    if debug
        println(res)
        println(mxh)
        println(mxh0)
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