using FUSE
using Test

@testset "nbi_power_weighted_fractions" begin
    t = collect(range(0.0, 5.0, 501))
    F = repeat([0.6, 0.3, 0.1], 1, length(t))
    P = fill(2e6, length(t))
    # beam off in the window: no fractions
    @test FUSE.nbi_power_weighted_fractions(t, zero(P), F, 1.0, 4.0) === nothing
    # only times with power count: the E/3 share of beam-off samples is ignored
    F2 = copy(F)
    F2[:, t.<2.0] .= [0.2, 0.2, 0.6]
    P2 = ifelse.(t .< 2.0, 0.0, 2e6)
    @test FUSE.nbi_power_weighted_fractions(t, P2, F2, 0.0, 5.0) ≈ [0.6, 0.3, 0.1]
    # weighted by the power launched: 1 MW of a then 3 MW of b for equal durations
    a, b = [0.7, 0.2, 0.1], [0.5, 0.3, 0.2]
    F3 = ifelse.(t' .< 2.5, a, b)
    P3 = ifelse.(t .< 2.5, 1e6, 3e6)
    @test FUSE.nbi_power_weighted_fractions(t, P3, F3, 0.0, 5.0) ≈ (a .+ 3 .* b) ./ 4 rtol = 1e-2
end
