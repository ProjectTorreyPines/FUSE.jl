using FUSE
using Test

@testset "nb_component_power_fractions" begin
    dd = IMAS.dd()
    dd.global_time = 1.0
    nbu = resize!(dd.nbi.unit, 1)[1]
    nbu.name = "30L"
    fbcur = [0.6, 0.3, 0.1]
    @test FUSE.nb_component_power_fractions(nbu, fbcur, :beam_current_fraction) === fbcur
    # no beam_power_fraction: the power carried by the current fractions, P_k ∝ f_k E/k
    p = fbcur ./ (1:3)
    @test FUSE.nb_component_power_fractions(nbu, fbcur, :beam_power_fraction) ≈ p ./ sum(p)
    # beam_power_fraction, normalized
    nbu.beam_power_fraction.time = [0.0]
    nbu.beam_power_fraction.data = reshape([0.8, 0.3, 0.1], :, 1)
    @test FUSE.nb_component_power_fractions(nbu, fbcur, :beam_power_fraction) ≈ [0.8, 0.3, 0.1] ./ 1.2
    # an invalid one falls back to the current fractions
    nbu.beam_power_fraction.data = reshape([0.8, -0.3, 0.1], :, 1)
    @test FUSE.nb_component_power_fractions(nbu, fbcur, :beam_power_fraction) ≈ p ./ sum(p)
end
