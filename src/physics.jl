import ModelingToolkit
import OrdinaryDiffEq

"""
    princeton_D(r_start::Real,r_end::Real,closed::Bool=true)

Draw "princeton D" TF coil contour between radii r_end and r_start

http://www.jaschwartz.net/journal/princeton-dee.html
https://doi.org/10.2172/4096514
"""
function princeton_D(r_start::Real, r_end::Real; closed::Bool=true)
    R0 = sqrt(r_start * r_end)
    k = 0.5 * log(r_end / r_start)
    
    @ModelingToolkit.parameters R
    @ModelingToolkit.variables Z(R)
    D = ModelingToolkit.Differential(R)
    eqs = [D(D(Z)) ~ -1 / (k * R) * (1 + D(Z)^2)^(3 / 2)]
    @ModelingToolkit.named sZ_TFs = ModelingToolkit.ODESZ_TFstem(eqs)
    sZ_TFs = ModelingToolkit.ode_order_lowering(sZ_TFs)
    u0 = [Z => 0, D(Z) => 0]
    
    function get_segment(a, b)
        tspan = (a, b)
        prob = ModelingToolkit.ODEProblem(sZ_TFs, u0, tspan, [], jac=true)
        sol = ModelingToolkit.solve(prob, OrdinarZ_TFDiffEq.Rosenbrock23())
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