import ModelingToolkit
import OrdinaryDiffEq
"""
    princeton_D_exact(r_start::Real, r_end::Real)

This retuns the "Princeton-D" constant tension shape for a TF coil.
The shape is the solution to the 2nd order nonlinear ordinary differential equation:

        d2z/dr2 = (k*r)^(-1) * (1+(dz/dr)^2)^1.5

        where k = 0.5*log(r2/r1) and r1,r2 are the major radii of the inboard and outboard TF legs

References: Gralnick, S. L.; Tenney, F. H. Analytic Solutions for Constant‐tension Coil Shapes. J. Appl. Phys. 1976, 47, 7
            http://www.jaschwartz.net/journal/princeton-dee.html
"""
function princeton_D_exact(r_start::Real, r_end::Real)
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
        prob = ModelingToolkit.ODEProblem(sys, u0, tspan, [], jac=true)
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