using FUSE
using Test

@testset "ActorSimpleNB.axis_regularization" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    ρ = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    function nb(rho_reg)
        act.ActorSimpleNB.axis_regularization = rho_reg
        FUSE.ActorSimpleNB(dd, act)
        s = findfirst(:nbi, dd.core_sources.source).profiles_1d[]
        totals = (s.electrons.power_inside[end], s.total_ion_power_inside[end], s.electrons.particles_inside[end], s.torque_tor_inside[end], s.current_parallel_inside[end])
        return totals, copy(s.electrons.energy), copy(s.j_parallel)
    end
    # second-order one-sided derivative at rho=0 (exactly zero for an even quadratic), normalized
    slope(q) = abs(-3q[1] + 4q[2] - q[3]) / (ρ[3] - ρ[1]) / maximum(abs, q)
    T0, qe0, j0 = nb(0.0)
    T1, qe1, j1 = nb(0.05)
    # same totals
    for (t0, t1) in zip(T0, T1)
        @test t1 ≈ t0 rtol = 1e-10
    end
    # untouched outside the regularization region
    out = ρ .>= 0.05
    @test qe1[out] == qe0[out]
    @test j1[out] == j0[out]
    # (nearly) zero slope on axis
    @test slope(qe1) < slope(qe0) / 20
    @test slope(j1) < slope(j0) / 20
end
