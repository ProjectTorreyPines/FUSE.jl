import ModelingToolkit
import OrdinaryDiffEq
import LazySets
import Interpolations
import PolygonOps

#= =============== =#
#  Shape functions  #
#= =============== =#

function error_layer_shape(shape_function_index)
    error("layer.shape=$(shape_function_index) is invalid. Valid options are:
 -2: Offset & convex-hull
 -1: Offset
  1: Priceton D  (shape_parameters = [])
  2: rectangle   (shape_parameters = [height])
  3: tripple-arc (shape_parameters = [height, small_radius, mid_radius, small_coverage, mid_coverage])
  4: miller      (shape_parameters = [elongation, triangularity])
  5: spline      (shape_parameters = [hfact, rz...)
101: Priceton D  (shape_parameters = [z_offset])
102: rectangle   (shape_parameters = [height, z_offset])
103: tripple-arc (shape_parameters = [height, small_radius, mid_radius, small_coverage, mid_coverage, z_offset])
104: miller      (shape_parameters = [elongation, triangularity, z_offset])
105: spline      (shape_parameters = [hfact, rz..., z_offset])
")
end

function init_shape_parameters(shape_function_index, r_obstruction, z_obstruction, r_start, r_end, target_clearance)
    height = maximum(z_obstruction) - minimum(z_obstruction) + target_clearance * 2.0
    z_offset = (maximum(z_obstruction) + minimum(z_obstruction)) / 2.0
    shape_parameters = nothing
    if shape_function_index in [-1, -2]
        return nothing
    else
        shape_index_mod = mod(shape_function_index, 100)
        if shape_index_mod == 1
            shape_parameters = Real[]
        elseif shape_index_mod == 2
            shape_parameters = [height]
        elseif shape_index_mod == 3
            shape_parameters = [log10(height), log10(1E-3), log10(1E-3), log10(45), log10(45)]
        elseif shape_index_mod == 4
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
            shape_parameters = [elongation, (triup + tridown) / 2.0]
        elseif shape_index_mod == 5
            n = 2
            R = range(r_start, r_end, length = 2 + n)[2:end-1]
            Z = range(height / 2.0, height / 2.0, length = 2 + n)[2:end-1]
            shape_parameters = Float64[0.8]
            for (r, z) in zip(R, Z)
                append!(shape_parameters, [r, z])
            end
        end
    end
    if shape_parameters === nothing
        error(shape_function_index)
    end
    if shape_function_index > 100
        push!(shape_parameters, z_offset)
    end
    return shape_parameters
end

function add_z_offset_parameter(func)
    function (args...)
        R, Z = func(args[1:end-1]...)
        Z .+= args[end]
        return R, Z
    end
end

function shape_function(shape_function_index)
    func = nothing
    if shape_function_index in [-1, -2]
        return nothing
    else
        shape_index_mod = mod(shape_function_index, 100)
        if shape_index_mod in 1
            func = princeton_D
        elseif shape_index_mod == 2
            func = (r_start, r_end, height) -> rectangle_shape(r_start, r_end, height; n_points = 100)
        elseif shape_index_mod == 3
            func = tripple_arc
        elseif shape_index_mod == 4
            func = miller_Rstart_Rend
        elseif shape_index_mod == 5
            func = spline_shape
        end
    end
    if func === nothing
        error_layer_shape(shape_function_index)
    end
    if shape_function_index > 100
        func = add_z_offset_parameter(func)
    end
    return func
end

"""
    optimize_shape(r_obstruction, z_obstruction, target_clearance, func, r_start, r_end, shape_parameters; verbose=false, time_limit=60)

Find shape parameters that generate smallest shape and target clearance from an obstruction
"""
function optimize_shape(r_obstruction, z_obstruction, target_clearance, func, r_start, r_end, shape_parameters; verbose = false, time_limit = 60)

    if length(shape_parameters) in [0, 1]
        func(r_start, r_end, shape_parameters...)
        return shape_parameters
    end

    function cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, obstruction_area, target_clearance, func, r_start, r_end, shape_parameters)
        R, Z = func(r_start, r_end, shape_parameters...)

        # disregard near r_start and r_end where optimizer has no control and shape is allowed to go over obstruction
        index = abs.((R .- r_start) .* (R .- r_end)) .> (target_clearance / 1000)
        R = R[index]
        Z = Z[index]

        # minimize area  O(1)
        area = sum(abs.(diff(R) .* (Z[1:end-1] .+ Z[2:end])))
        cost_area = (area - obstruction_area) / obstruction_area

        # no polygon crossings  O(N)
        inpoly = [PolygonOps.inpolygon((r, z), rz_obstruction) for (r, z) in zip(R, Z)]
        cost_inside = sum(inpoly) / cost_area

        # target clearance  O(1)
        minimum_distance = minimum_distance_two_shapes(R, Z, r_obstruction, z_obstruction)
        cost_min_clearance = (minimum_distance - target_clearance) / target_clearance

        # favor up/down symmetric solutions
        cost_up_down_symmetry = abs(sum(Z) / sum(abs.(Z)))

        # favor smoothness and avoid convex shapes
        dR = diff(R)
        dZ = diff(Z)
        dRa = dR[1:end-1]
        dRb = dR[2:end]
        dZa = dZ[1:end-1]
        dZb = dZ[2:end]
        sin_theta = (dRa .* dZb .- dZa .* dRb) ./ (sqrt.(dRa .^ 2.0 .+ dZa .^ 2.0) .* sqrt.(dRb .^ 2.0 + dZb .^ 2.0))
        cost_concavity = 0
        if any(sin_theta .> 0)
            cost_concavity = sum(abs.(sin_theta[sin_theta.>0])) / length(sin_theta)
        end
        cost_smoothness = sum(abs.(sin_theta)) / length(sin_theta)

        # return cost
        return cost_min_clearance^2 + 0.5 * cost_area^2 + cost_inside^2 + 1E-3 * cost_up_down_symmetry^2 + 10 * cost_smoothness^2 + 100 * cost_concavity^2
    end

    rz_obstruction = collect(zip(r_obstruction, z_obstruction))
    obstruction_area = sum(abs.(diff(r_obstruction) .* (z_obstruction[1:end-1] .+ z_obstruction[2:end])))
    initial_guess = copy(shape_parameters)
    # res = optimize(shape_parameters-> cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, obstruction_area, target_clearance, func, r_start, r_end, shape_parameters),
    #                initial_guess, Newton(), Optim.Options(time_limit=time_limit); autodiff=:forward)
    res = Optim.optimize(shape_parameters -> cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, obstruction_area, target_clearance, func, r_start, r_end, shape_parameters),
        initial_guess, length(shape_parameters) == 1 ? Optim.BFGS() : Optim.NelderMead(), Optim.Options(time_limit = time_limit); autodiff = :forward)
    if verbose
        println(res)
    end
    shape_parameters = Optim.minimizer(res)
    R, Z = func(r_start, r_end, shape_parameters...)
    # plot(func(r_start, r_end, initial_guess...);markershape=:+)
    # plot!(r_obstruction,z_obstruction)
    # display(plot!(R,Z;markershape=:x,aspect_ratio=:equal))
    return shape_parameters
end

"""
    princeton_D(r_start::Real, r_end::Real)

"Princeton D" contour between radii r_end and r_start

http://www.jaschwartz.net/journal/princeton-dee.html
https://doi.org/10.2172/4096514

layer[:].shape = 1
"""
function princeton_D(r_start::Real, r_end::Real)
    R0 = sqrt(r_start * r_end)
    k = 0.5 * log(r_end / r_start)

    ModelingToolkit.@parameters R
    ModelingToolkit.@variables Z(R)
    D = ModelingToolkit.Differential(R)
    eqs = [D(D(Z)) ~ -1 / (k * R) * (1 + D(Z)^2)^(3 / 2)]
    ModelingToolkit.@named sys = ModelingToolkit.ODESystem(eqs)
    sys = ModelingToolkit.ode_order_lowering(sys)
    u0 = [Z => 0, D(Z) => 0]

    function get_segment(a, b)
        tspan = (a, b)
        prob = ModelingToolkit.ODEProblem(sys, u0, tspan, [], jac = true)
        sol = ModelingToolkit.solve(prob, OrdinaryDiffEq.Rosenbrock23())
        return sol.t, [u[2] for u in sol.u]
    end

    segment1 = get_segment(R0, r_start)
    segment2 = get_segment(R0, r_end)

    Z02 = segment2[2][end]
    segment1[2] .-= Z02
    segment2[2] .-= Z02

    R = vcat(reverse(segment1[1])[1:end-1], segment2[1][1:end-1], reverse(segment2[1])[1:end-1], segment1[1], segment1[1][end])
    Z = vcat(reverse(segment1[2])[1:end-1], segment2[2][1:end-1], -reverse(segment2[2])[1:end-1], -segment1[2], segment1[2][end])

    return R, Z
end

"""
    rectangle_shape(r_start::Real, r_end::Real, z_low::Real, z_high::Real; n_points = 5)

Asymmetric rectangular contour
layer[:].shape = 2
"""
function rectangle_shape(r_start::Real, r_end::Real, z_low::Real, z_high::Real; n_points = 5)
    if n_points == 5
        R = [r_start, r_end, r_end, r_start, r_start]
        Z = [z_low, z_low, z_high, z_high, z_low]
    else
        R = vcat(range(r_start, r_end, length = n_points), range(r_end, r_end, length = n_points)[2:end], range(r_end, r_start, length = n_points)[2:end], range(r_start, r_start, length = n_points)[2:end], r_start)
        Z = vcat(range(z_low, z_low, length = n_points), range(z_low, z_high, length = n_points)[2:end], range(z_high, z_high, length = n_points)[2:end], range(z_high, z_low, length = n_points)[2:end], z_low)
    end
    return R, Z
end

"""
    rectangle_shape(r_start::Real, r_end::Real, height::Real; n_points = 5)

Symmetric rectangular contour
"""
function rectangle_shape(r_start::Real, r_end::Real, height::Real; n_points = 5)
    Δ = height / 2.0
    return rectangle_shape(r_start, r_end, -Δ, Δ; n_points)
end

"""
    function tripple_arc(r_start::Real,
        r_end::Real,
        height::Real,
        small_radius::Real,
        mid_radius::Real,
        small_coverage::Real,
        mid_coverage::Real;
        min_small_radius_fraction::Real=0.5,
        min_mid_radius_fraction::Real=min_small_radius_fraction*2.0,
        n_points::Int=400)

TrippleArc contour
Angles are in degrees
height, small_radius, mid_radius, small_coverage, mid_coverage are 10^exponent (to ensure positiveness)
"""
function tripple_arc(r_start::Real,
    r_end::Real,
    height::Real,
    small_radius::Real,
    mid_radius::Real,
    small_coverage::Real,
    mid_coverage::Real;
    min_small_radius_fraction::Real = 0.25,
    min_mid_radius_fraction::Real = min_small_radius_fraction * 2.0,
    n_points::Int = 400)

    height = 10^height / 2.0
    small_radius = 10^small_radius + height * min_small_radius_fraction
    mid_radius = 10^mid_radius + height * min_mid_radius_fraction
    small_coverage = 10^small_coverage * pi / 180
    mid_coverage = 10^mid_coverage * pi / 180

    asum = small_coverage + mid_coverage
    n_points = floor(Int, n_points / 4)

    # small arc
    theta = LinRange(0, small_coverage, n_points)
    small_arc_R = r_start .+ small_radius .* (1 .- cos.(theta))
    small_arc_Z = height .+ small_radius .* sin.(theta)

    # mid arc
    theta = LinRange(small_coverage, asum, n_points)
    mid_arc_R = small_arc_R[end] .+ mid_radius .* (cos.(small_coverage) .- cos.(theta))
    mid_arc_Z = small_arc_Z[end] .+ mid_radius .* (sin.(theta) .- sin.(small_coverage))

    # large arc
    theta = LinRange(theta[end], pi, n_points)
    large_radius = mid_arc_Z[end] / sin(pi - asum)
    large_arc_R = mid_arc_R[end] .+ large_radius .* (cos.(pi .- theta) .- cos.(pi .- asum))
    large_arc_Z = mid_arc_Z[end] .- large_radius .* (sin(asum) .- sin.(pi .- theta))

    R = vcat(small_arc_R, mid_arc_R[2:end], large_arc_R[2:end])
    R = vcat(R, reverse(R)[2:end])
    Z = vcat(small_arc_Z, mid_arc_Z[2:end], large_arc_Z[2:end])
    Z = vcat(Z, -reverse(Z)[2:end])

    # Add vertical
    R = vcat(LinRange(r_start, r_start, n_points), R)
    Z = vcat(LinRange(-height, height, n_points), Z)

    # Resize to ensure r_start to r_end
    factor = (r_end - r_start) / (maximum(R) - minimum(R))
    Z = Z .* factor
    R = (R .- minimum(R)) .* factor .+ r_start

    return R, Z
end

"""
    miller(R0, inverse_aspect_ratio, elongation, triangularity, n_points)

Miller contour
layer[:].shape = 4
"""
function miller(R0::Real, rmin_over_R0::Real, elongation::Real, triangularity::Real; n_points::Int = 401)
    θ = range(0, 2 * pi, length = n_points)
    # bound triangularity
    while abs(triangularity) > 1.0
        if triangularity < 1.0
            triangularity = -2.0 - triangularity
        else
            triangularity = 2.0 - triangularity
        end
    end
    δ₀ = asin(triangularity)
    R = R0 * (1 .+ rmin_over_R0 .* cos.(θ .+ δ₀ * sin.(θ)))
    Z = R0 * (rmin_over_R0 * elongation * sin.(θ))
    R[end] = R[1]
    Z[end] = Z[1]
    return R, Z
end

"""
    miller_Rstart_Rend(r_start, r_end, elongation, triangularity, n_points)

Miller contour
"""
function miller_Rstart_Rend(r_start::Real, r_end::Real, elongation::Real, triangularity::Real; n_points::Int = 401)
    return miller((r_end + r_start) / 2.0, (r_end - r_start) / (r_end + r_start), elongation, triangularity; n_points)
end

"""
    spline_shape(r::Real, z::Real; n_points::Int=101)

Spline contour
"""
function spline_shape(r::Vector{T}, z::Vector{T}; n_points::Int = 101) where {T<:Real}
    r = vcat(r[1], r[1], r, r[end], r[end])
    z = vcat(0, z[1] / 2, z, z[end] / 2, 0)
    d = cumsum(sqrt.(vcat(0, diff(r)) .^ 2.0 .+ vcat(0, diff(z)) .^ 2.0))

    itp_r = Interpolations.interpolate(d, r, Interpolations.FritschButlandMonotonicInterpolation())
    itp_z = Interpolations.interpolate(d, z, Interpolations.FritschButlandMonotonicInterpolation())

    D = LinRange(d[1], d[end], n_points)
    R, Z = itp_r.(D), itp_z.(D)
    R[end] = R[1]
    Z[end] = Z[1]
    return R, Z
end

function spline_shape(r_start::Real, r_end::Real, hfact::Real, rz...; n_points = 101)
    rz = collect(rz)
    R = rz[1:2:end]
    Z = rz[2:2:end]
    hfact_max = 0.5 + (minimum(R) - r_start) / (r_end - r_start) / 2.0
    hfact = min(abs(hfact), hfact_max)
    h = maximum(Z) * hfact
    r = vcat(r_start, R, r_end, reverse(R), r_start)
    z = vcat(h, Z, 0, -reverse(Z), -h)
    return spline_shape(r, z; n_points = n_points)
end


function xy_polygon(x, y)
    if (x[1] != x[end]) && (x[1] ≈ x[end])
        x[end] = x[1]
    end
    if (y[1] != y[end]) && (y[1] ≈ y[end])
        y[end] = y[1]
    end
    if (x[1] != x[end]) || (y[1] != y[end])
        push!(x, x[1])
        push!(y, y[1])
    end
    coords = [collect(map(collect, zip(x, y)))]
    return LibGEOS.Polygon(coords)
end