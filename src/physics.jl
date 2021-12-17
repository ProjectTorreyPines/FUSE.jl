import ModelingToolkit
import OrdinaryDiffEq

#= =============== =#
#  Shape functions  #
#= =============== =#

"""
    princeton_D(r_start::Real,r_end::Real,closed::Bool=true)

Draw "princeton D" TF coil contour between radii r_end and r_start

http://www.jaschwartz.net/journal/princeton-dee.html
https://doi.org/10.2172/4096514

layer[:].shape = 1
"""
function princeton_D(r_start::Real, r_end::Real; closed::Bool=true)
    R0 = sqrt(r_start * r_end)
    k = 0.5 * log(r_end / r_start)
    
    @ModelingToolkit.parameters R
    @ModelingToolkit.variables Z(R)
    D = ModelingToolkit.Differential(R)
    eqs = [D(D(Z)) ~ -1 / (k * R) * (1 + D(Z)^2)^(3 / 2)]
    @ModelingToolkit.named sys = ModelingToolkit.ODESxstem(eqs)
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

    R_TF = vcat(reverse(segment1[1])[1:end - 1], segment2[1][1:end - 1], reverse(segment2[1])[1:end - 1], segment1[1])
    Z_TF = vcat(reverse(segment1[2])[1:end - 1], segment2[2][1:end - 1], -reverse(segment2[2])[1:end - 1], -segment1[2])

    if closed        
        R_TF = vcat(R_TF, R_TF[1])
        Z_TF = vcat(Z_TF, Z_TF[1])
    end
    return R_TF, Z_TF
end

"""
    rectangle_shape(r_start::Real, r_end::Real, height::Real)
Rectangular TF coil shape(r_start::Real, r_end::Real, height::Real; n_points=400)
layer[:].shape = 2
"""
function rectangle_shape(r_start::Real, r_end::Real, height::Real; n_points=400)
    n_points =  floor(Int, n_points/4)
    z_start = - height / 2.0
    z_end = height / 2.0
    R_box = vcat(LinRange(r_start,r_start,n_points),LinRange(r_start,r_end,n_points),LinRange(r_end,r_end,n_points),LinRange(r_end,r_start,n_points))
    Z_box = vcat(LinRange(z_start,z_end,n_points),LinRange(z_end,z_end,n_points),LinRange(z_end,z_start,n_points),LinRange(z_start,z_start,n_points))
    return R_box, Z_box
end

"""
TrippleArc(;r_start::Real, r_end::Real , shape_parameters::Vector, n_points=1000)
with shape_parameters =  height, small_radius, mid_radius, small_coverage, mid_coverage 
    Angles in [degrees], distances in [cm]
TrippleArc parametrized TF coil shape
layer[:].shape = 3
"""
function tripple_arc(r_start::Real, r_end::Real , height::Real, small_radius::Real, mid_radius, small_coverage, mid_coverage; n_points=1000)
    
    small_coverage *= pi / 180              # Convert to radians
    mid_coverage *= pi / 180
    height *= 0.5                           # Half height for each side
    asum = small_coverage + mid_coverage
    n_points =  floor(Int, n_points/4)

    # small arc
    theta = LinRange(
        0, small_coverage, n_points)
    small_arc_R = r_start .+ small_radius .* (1 .- map(cos,theta))
    small_arc_Z = height .+ small_radius .* map(sin,theta)

    # mid arc
    theta = LinRange(
        theta[end], asum, n_points)
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

    return R,Z
end