import ModelingToolkit
import OrdinaryDiffEq
import LazySets

#= =============== =#
#  Shape functions  #
#= =============== =#

function init_shape_parameters(shape_function_index, r_obstruction, z_obstruction, r_start, r_end, target_clearance)
    height = maximum(z_obstruction) - minimum(z_obstruction) + target_clearance * 2
    if shape_function_index in [1, -1, -2]
        shape_parameters = Real[]
    elseif shape_function_index == 2
        shape_parameters = [height]
    elseif shape_function_index == 3
        shape_parameters = [height, height*0.1, height*0.25, 45, 90]
    elseif shape_function_index == 4
        shape_parameters = [height/(r_end-r_start), 0.0]
    end
    return shape_parameters
end

function shape_function(shape_index)
    if shape_index in [-1, -2]
        return nothing
    elseif shape_index  == 1
        return princeton_D
    elseif shape_index == 2
        return (r_start, r_end, height) -> rectangle_shape(r_start, r_end, height; n_points = 100)
    elseif shape_index == 3
        return  tripple_arc
    elseif shape_index == 4
        return miller_Rstart_Rend
    else
        error("layer.shape=$(shape) is invalid. Valid options are:
1: Priceton D  (shape_parameters = [])
2: rectangle   (shape_parameters = [height])
3: tripple-arc (shape_parameters = [height, small_radius, mid_radius, small_coverage, mid_coverage])
4: miller      (shape_parameters = [elongation, triangularity])
")
    end
end

"""
    optimize_shape_clearance(r_obstruction, z_obstruction, target_clearance, func, r_start, r_end, shape_parameters; verbose=false, time_limit=60)

Find shape parameters that generate smallest shape and target clearance from an obstruction
"""
function optimize_shape(r_obstruction, z_obstruction, target_clearance, func, r_start, r_end, shape_parameters; verbose=false, time_limit=60)

    if length(shape_parameters) == 0
        func(r_start, r_end, shape_parameters...)
        return shape_parameters 
    end

    function cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, obstruction_area, target_clearance, func, r_start, r_end, shape_parameters)
        R, Z = func(r_start, r_end, shape_parameters...)
        
        # disregard near r_start and r_end where optimizer has no control and shape is allowed to go over obstruction
        index = abs.((R .- r_start) .* (R .- r_end)) .> (target_clearance/1000)
        R = R[index]
        Z = Z[index]

        # no polygon crossings!  O(N)
        inpoly = [PolygonOps.inpolygon((r, z), rz_obstruction) for (r,z) in zip(R, Z)]
        cost_inside = sum(inpoly)

        # minimize area  O(1)
        coil_area = sum(abs.(diff(R) .* (Z[1:end-1] .+ Z[2:end])))
        cost_area = (coil_area - obstruction_area) / obstruction_area
        
        # target clearance  O(1)
        minimum_distance = minimum_distance_two_shapes(R, Z, r_obstruction, z_obstruction)
        cost_min_clearance = (minimum_distance - target_clearance) / target_clearance
        
        # return cost
        return cost_min_clearance^2 + 1E-1 * cost_area^2 + cost_inside^2
    end

    rz_obstruction = collect(zip(r_obstruction, z_obstruction))
    obstruction_area =  sum(abs.(diff(r_obstruction) .* (z_obstruction[1:end-1] .+ z_obstruction[2:end]) ))
    initial_guess = copy(shape_parameters)
    res = optimize(shape_parameters-> cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, obstruction_area, target_clearance, func, r_start, r_end, shape_parameters),
                   initial_guess, length(shape_parameters)==1 ? BFGS() : NelderMead(), Optim.Options(time_limit=time_limit); autodiff=:forward)
    if verbose
        println(res)
    end
    shape_parameters = Optim.minimizer(res)
    return shape_parameters
end

"""
    miller(R0, inverse_aspect_ratio, elongation, triangularity, n_points)

Miller contour
"""
function miller(R0, rmin_over_R0, elongation, triangularity; n_points = 401)
    θ = range(0, 2*pi, length=n_points)
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
function miller_Rstart_Rend(r_start, r_end, elongation, triangularity; n_points = 401)
    return miller((r_end + r_start) / 2.0, (r_end - r_start) / (r_end + r_start), elongation, triangularity; n_points)
end

"""
    princeton_D(r_start::Real,r_end::Real,closed::Bool=true)

"Princeton D" contour between radii r_end and r_start

http://www.jaschwartz.net/journal/princeton-dee.html
https://doi.org/10.2172/4096514

layer[:].shape = 1
"""
function princeton_D(r_start::Real, r_end::Real)
    R0 = sqrt(r_start * r_end)
    k = 0.5 * log(r_end / r_start)
    
    @ModelingToolkit.parameters R
    @ModelingToolkit.variables Z(R)
    D = ModelingToolkit.Differential(R)
    eqs = [D(D(Z)) ~ -1 / (k * R) * (1 + D(Z)^2)^(3 / 2)]
    @ModelingToolkit.named sys = ModelingToolkit.ODESystem(eqs)
    sys = ModelingToolkit.ode_order_lowering(sys)
    u0 = [Z => 0, D(Z) => 0]
    
    function get_segment(a, b)
        tspan = (a, b)
        prob = ModelingToolkit.ODEProblem(sys, u0, tspan, [], jac=true)
        sol = ModelingToolkit.solve(prob, OrdinaryDiffEq.Rosenbrock23())
        return sol.t, [u[2] for u in sol.u]
    end

    segment1 = get_segment(R0, r_start)
    segment2 = get_segment(R0, r_end)

    Z02 = segment2[2][end]
    segment1[2] .-= Z02
    segment2[2] .-= Z02

    R = vcat(reverse(segment1[1])[1:end - 1], segment2[1][1:end - 1], reverse(segment2[1])[1:end - 1], segment1[1], segment1[1][end])
    Z = vcat(reverse(segment1[2])[1:end - 1], segment2[2][1:end - 1], -reverse(segment2[2])[1:end - 1], -segment1[2], segment1[2][end])

    return R, Z
end


"""
    rectangle_shape(r_start::Real, r_end::Real, z_low::Real, z_high::Real)

Asymmetric rectangular contour
"""
function rectangle_shape(r_start::Real, r_end::Real, z_low::Real, z_high::Real; n_points = 5)
    if n_points == 5
        R = [r_start, r_end, r_end, r_start, r_start]
        Z = [z_low, z_low, z_high, z_high, z_low]
    else
        R = vcat(range(r_start, r_end, length=n_points),range(r_end, r_end, length=n_points)[2:end],range(r_end, r_start, length=n_points)[2:end],range(r_start, r_start, length=n_points)[2:end],r_start)
        Z = vcat(range(z_low, z_low, length=n_points),range(z_low, z_high, length=n_points)[2:end],range(z_high, z_high, length=n_points)[2:end],range(z_high, z_low, length=n_points)[2:end],z_low)
    end
    return R, Z
end

"""
    rectangle_shape(r_start::Real, r_end::Real, height::Real)

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
                         n_points::Int=1000)

TrippleArc contour
Angles are in degrees
"""
function tripple_arc(r_start::Real,
                     r_end::Real,
                     height::Real,
                     small_radius::Real,
                     mid_radius::Real,
                     small_coverage::Real,
                     mid_coverage::Real;
                     min_small_radius_fraction::Real=0.1,
                     min_mid_radius_fraction::Real=0.25,
                     n_points::Int=1000)

    height = abs(height) * 0.5
    small_radius = abs(small_radius) + height * min_small_radius_fraction
    mid_radius = abs(mid_radius) + height * min_mid_radius_fraction
    small_coverage = abs(small_coverage) * pi / 180
    mid_coverage =  abs(mid_coverage) * pi / 180

    asum = small_coverage + mid_coverage
    n_points =  floor(Int, n_points/4)

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
    R = vcat(R,reverse(R)[2:end])
    Z = vcat(small_arc_Z, mid_arc_Z[2:end], large_arc_Z[2:end])
    Z = vcat(Z, -reverse(Z)[2:end])
    
    # Add vertical
    R = vcat(LinRange(r_start, r_start, n_points), R)
    Z = vcat(LinRange(-height, height ,n_points), Z)
    
    # Resize to ensure r_start to r_end
    factor = (r_end - r_start) / (maximum(R) - minimum(R))
    Z = Z .* factor
    R = (R .- minimum(R)) .* factor .+ r_start

    return R, Z
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