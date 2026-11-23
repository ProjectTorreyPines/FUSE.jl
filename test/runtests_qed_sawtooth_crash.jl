using FUSE
using Test

@testset "sawtooth_reset_q" begin
    ρ = collect(range(0, 1, 101))
    q_s, dq = 1.025, 0.03

    # nothing below q_s: no region, q unchanged
    q = 1.1 .+ 2 .* ρ .^ 2
    r = FUSE.sawtooth_reset_q(ρ, q, q_s, dq; ramp=0.01)
    @test isempty(r.regions) && r.q_new == q

    # one off-axis dip: reset to q_s at its minimum, q unchanged outside [ρ_in, ρ_out] and on the axis
    q = 1.04 .- 0.05 .* exp.(-((ρ .- 0.3) ./ 0.08) .^ 2) .+ 2 .* ρ .^ 4
    q[1] = q[2]
    r = FUSE.sawtooth_reset_q(ρ, q, q_s, dq; ramp=0.0)
    @test length(r.regions) == 1
    g = r.regions[1]
    out = (ρ .<= g.ρ_in) .| (ρ .>= g.ρ_out)
    @test r.q_new[out] == q[out]
    # q_s is set at the dip centre, refined between nodes: the node minimum is just above it
    @test q_s <= minimum(r.q_new) < q_s + 1e-4

    # ramp: a dip just below q_s is only partly reset, continuously
    q2 = q .+ (q_s - 0.002 - minimum(q))
    r2 = FUSE.sawtooth_reset_q(ρ, q2, q_s, dq; ramp=0.01)
    @test r2.regions[1].w ≈ 0.2
    @test minimum(r2.q_new) < q_s

    # multi-dip profiles: the reset never takes q past q_s + dq where it was below it
    for q0 in 0.97:0.01:1.05, k in 1:4
        q = q0 .+ 2 .* ρ .^ 2 .+ 0.02 .* sin.(2π * k .* ρ)
        q[1] = q[2]
        r = FUSE.sawtooth_reset_q(ρ, q, q_s, dq; ramp=0.01)
        @test all(r.q_new[q.<q_s+dq] .<= q_s + dq + 1e-12)
    end
end
