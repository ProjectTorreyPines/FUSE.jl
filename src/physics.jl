import ModelingToolkit
import OrdinaryDiffEq

"""
    princeton_D(R1::Real,R2::Real,closed::Bool=true)

Draw "princeton D" TF coil contour between radii R2 and R1

http://www.jaschwartz.net/journal/princeton-dee.html
https://doi.org/10.2172/4096514
"""
function princeton_D(R1::Real, R2::Real; closed::Bool=true)
    R0 = sqrt(R1 * R2)
    k = 0.5 * log(R2 / R1)
    
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

    segment1 = get_segment(R0, R1)
    segment2 = get_segment(R0, R2)

    Z02 = segment2[2][end]
    segment1[2] .-= Z02
    segment2[2] .-= Z02

    x = vcat(reverse(segment1[1])[1:end - 1], segment2[1][1:end - 1], reverse(segment2[1])[1:end - 1], segment1[1])
    y = vcat(reverse(segment1[2])[1:end - 1], segment2[2][1:end - 1], -reverse(segment2[2])[1:end - 1], -segment1[2])

    if closed        
        x = vcat(x, x[1])
        y = vcat(y, y[1])
    end
    return x, y
end