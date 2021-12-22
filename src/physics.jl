import ModelingToolkit
import OrdinaryDiffEq

#= =============== =#
#  Shape functions  #
#= =============== =#

"""
    miller(R0, inverse_aspect_ratio, elongation, triangularity, n_points)

Miller contour
"""
function miller(R0, rmin_over_R0, elongation, triangularity, n_points)
    θ = range(0, 2pi, length=n_points)
    δ₀ = asin(triangularity)
    R = R0 * (1 .+ rmin_over_R0 .* cos.(θ .+ δ₀ * sin.(θ)))
    Z = R0 * (rmin_over_R0 * elongation * sin.(θ))
    return R, Z
end

"""
    miller_Rstart_Rend(r_start, r_end, elongation, triangularity, n_points)

Miller contour
"""
function miller_Rstart_Rend(r_start, r_end, elongation, triangularity, n_points)
    return miller((r_end + r_start) / 2.0, (r_end - r_start) / (r_end + r_start), elongation, triangularity, n_points)
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

    R = vcat(reverse(segment1[1])[1:end - 1], segment2[1][1:end - 1], reverse(segment2[1])[1:end - 1], segment1[1])
    Z = vcat(reverse(segment1[2])[1:end - 1], segment2[2][1:end - 1], -reverse(segment2[2])[1:end - 1], -segment1[2])

    return R, Z
end

"""
    rectangle_shape(r_start::Real, r_end::Real, z_low::Real, z_high::Real)

Asymmetric rectangular contour
"""
function rectangle_shape(r_start::Real, r_end::Real, z_low::Real, z_high::Real)
    return [r_start, r_end, r_end, r_start, r_start], [z_low, z_low, z_high, z_high, z_low]
end

"""
    rectangle_shape(r_start::Real, r_end::Real, height::Real)

Symmetric rectangular contour
"""
function rectangle_shape(r_start::Real, r_end::Real, height::Real)
    Δ = height / 2.0
    return rectangle_shape(r_start, r_end, -Δ, Δ)
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
                     n_points::Int=1000)
    
    small_coverage *= pi / 180              # Convert to radians
    mid_coverage *= pi / 180
    height *= 0.5                           # Half height for each side
    asum = small_coverage + mid_coverage
    n_points =  floor(Int, n_points/4)

    # small arc
    theta = LinRange(0, small_coverage, n_points)
    small_arc_R = r_start .+ small_radius .* (1 .- map(cos,theta))
    small_arc_Z = height .+ small_radius .* map(sin,theta)

    # mid arc
    theta = LinRange(theta[end], asum, n_points)
    mid_arc_R = small_arc_R[end] .+ mid_radius .* 
        (map(cos,small_coverage) .- map(cos,theta))
        mid_arc_Z = small_arc_Z[end] .+ mid_radius * 
        (map(sin,theta) .- map(sin,small_coverage))

    # large arc
    large_radius = (mid_arc_Z[end]) / sin(pi - asum)
    theta = LinRange(theta[end], pi, n_points)
    large_arc_R = mid_arc_R[end] .+ large_radius .* 
        (map(cos,pi .- theta) .- map(cos,pi .- asum))
        large_arc_Z = mid_arc_Z[end] .- large_radius .* 
        (sin(asum) .- map(sin, pi .- theta))

    R = vcat(small_arc_R, mid_arc_R[2:end], large_arc_R[2:end])
    R = vcat(R,reverse(R)[2:end])
    Z = vcat(small_arc_Z, mid_arc_Z[2:end], large_arc_Z[2:end])
    Z = vcat(Z, -reverse(Z)[2:end])
    
    # Add vertical
    R = vcat(LinRange(r_start, r_start, n_points), R)
    Z = vcat(LinRange(-height, height ,n_points), Z)
    
    # Resize to ensure r_start to r_end
    Z = Z./(maximum(R)-minimum(R)).*(r_end - r_start)
    R = (R .- minimum(R)) ./ (maximum(R) - minimum(R)) .* (r_end - r_start) .+ r_start

    return R, Z
end