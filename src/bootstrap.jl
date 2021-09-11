import Quadrature

abstract type AbstractBootstrapModel end

function bootstrap_fraction(bsm::AbstractBootstrapModel, betap, epsilon, St, Sn, Sj, Zeff)
    return betaP * collisionless_bootstrap(bsm, epsilon, St, Sn, Sj, Zeff) / sqrt(epsilon)
end

function bootstrap_fraction(fraction::Real, args...)
    return fraction
end

"""
Simple model for collisionless bootstrap coefficient predicted from
Pomphrey, N. Bootstrap dependence on plasma profile parameters. PPPL, 1992.
"""
struct CBSPomphrey <: AbstractBootstrapModel
end

function collisionless_bootstrap(bsm::CBSPomphrey, epsilon, St, Sn, Sj, Zeff)
    
    Sp = St + Sn
    
    Zp1 = (13.93 + 7.89 * Zeff) / (5.55 * Zeff - 1.0)
    Zt1 = -(34.37 + 0.63 * Zeff) / (8.11 * Zeff - 1.0)
    Zp2 = -(5.11 - 1.17 * Zeff) / (2.21 * Zeff - 1.0)
    Zt2 = (12.51 - 0.34 * Zeff) / (2.68 * Zeff - 1.0)
    
    g1_integrand = Quadrature.QuadratureProblem((x, p) -> x^1.25 * (1.0 - x)^(Sp - 1.0) / (1.0 - (1.0 - x)^(Sj + 1.0)), 0, 1)
    
    g2_integrand = Quadrature.QuadratureProblem((x, p) -> x^1.50 * (1.0 - x)^(Sp - 1.0) / (1.0 - (1.0 - x)^(Sj + 1.0)), 0, 1)
    
    I_g1 = solve(g1_integrand, Quadrature.QuadGKJL())
    I_g2 = solve(g2_integrand, Quadrature.QuadGKJL())
    
    G1 = 0.25 * (Sp + 1.0) * I_g1.u
    G2 = 0.25 * (Sp + 1.0) * I_g2.u
    
    cbs = (Sp * Zp1 + St * Zt1) * G1 + (Sp * Zp2 + St * Zt2) * G2 * sqrt(epsilon)
    
    return cbs
end

"""
Simple model for collisionless bootstrap coefficient predicted from
Gi et al., Fus. Eng. Design 89 2709 (2014)
"""
struct CBSGi <: AbstractBootstrapModel
end

function collisionless_bootstrap(bsm::CBSGi, epsilon, St, Sn, Sj, Zeff)

  Sp = St + Sn

  Cbs = 0.382 * epsilon^(-0.242) * Sp^0.974 * St^(-0.416) * Zeff^0.178

  return Cbs
end

"""
Simple model for collisionless bootstrap coefficient predicted from
Wilson et al., Nucl. Fusion 32 257 (1992)
"""
struct CBSWilson <: AbstractBootstrapModel
end

function collisionless_bootstrap(bsm::CBSWilson, epsilon, St, Sn, Sj, Zeff)

  Sp = Sn + St

  sqrt_eps = sqrt(epsilon)

  b = [1., Sp, St, Sp * St, sqrt_eps,
       Sp * sqrt_eps, St * sqrt_eps,
       Sp * St * sqrt_eps, epsilon,
       Sp * epsilon, St * epsilon,
       Sp * St * epsilon]

  a = [1.41 * (1 - 0.28 * sqrt(Sj)) * (1 + 0.12 / Zeff),
       0.36 * (1 - 0.59 * sqrt(Sj)) * (1 + 0.8 / Zeff),
       -0.27 * (1 - 0.47 * sqrt(Sj)) * (1 + 3. / Zeff),
       0.0053 * (1 + 5. / Zeff),
       -0.93 * (1. - 0.34 * sqrt(Sj)) * (1 + 0.15 / Zeff),
       -0.26 * (1. - 0.57 * sqrt(Sj)) * (1 - 0.27 * Zeff),
       0.064 * (1. - 0.6 * Sj + 0.15 * Sj^2.) * (1 + 7.6 / Zeff),
       -0.0011 * (1 + 9. / Zeff),
       -0.33 * (1. - Sj + 0.33 * Sj^2.),
       -0.26 * (1. - 0.87 / sqrt(Sj) - 0.16 * Sj),
       # *Corrected* PrevL ...-0.45alphaj
       -0.14 * (1. - 1.14 / sqrt(Sj) - 0.45 * sqrt(Sj)),
       -0.0069]


  Cbs = sum(a .* b)

  return Cbs

end
