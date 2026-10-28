using FUSE
using Test

@testset "ActorSimpleNB.axis_mirror" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    ρ = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    function nb(axis_mirror)
        act.ActorSimpleNB.axis_mirror = axis_mirror
        FUSE.ActorSimpleNB(dd, act)
        s = findfirst(:nbi, dd.core_sources.source).profiles_1d[]
        return s.electrons.power_inside[end] + s.total_ion_power_inside[end], s.electrons.energy
    end
    P0, qe0 = nb(false)
    P1, qe1 = nb(true)
    # same power, deposited with (nearly) zero slope on axis
    @test P1 ≈ P0 rtol = 1e-12
    slope(q) = abs(q[2] - q[1]) / (ρ[2] - ρ[1]) / maximum(q)
    @test slope(qe1) < slope(qe0) / 20
end
