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

@testset "regularize_on_axis! is linear" begin
    # components regularized separately must still sum to the regularized total
    ρ = range(0.0, 1.0, 101)
    vol = ρ .^ 2
    f = exp.(-(ρ .- 0.05) .^ 2 ./ 0.02)
    g = 1.0 .+ sqrt.(ρ)
    total = FUSE.regularize_on_axis!(f .+ 2g, ρ, vol, 0.1)
    @test FUSE.regularize_on_axis!(copy(f), ρ, vol, 0.1) .+ 2 .* FUSE.regularize_on_axis!(copy(g), ρ, vol, 0.1) ≈ total rtol = 1e-12
    # works in place on a view of a (group, component, rho) array, as used in ActorSimpleNB
    arr = zeros(2, 3, length(ρ))
    arr[2, 1, :] .= f
    FUSE.regularize_on_axis!(@view(arr[2, 1, :]), ρ, vol, 0.1)
    @test arr[2, 1, :] ≈ FUSE.regularize_on_axis!(copy(f), ρ, vol, 0.1)
end
