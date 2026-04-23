using FUSE
using Test

@testset "forward_ad flux matcher" begin
    # Initialize ITER case (lightweight scalar init)
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    # Configure for AD-compatible flux matching: TJLF + Hirshman-Sigmar
    act.ActorTGLF.model = :TJLF
    act.ActorNeoclassical.model = :hirshmansigmar
    act.ActorCoreTransport.model = :FluxMatcher
    act.ActorFluxMatcher.jacobian_method = :forward_ad
    act.ActorFluxMatcher.max_iterations = 3
    act.ActorFluxMatcher.evolve_pedestal = false

    # Save initial core_profiles for comparison
    cp1d_before = deepcopy(dd.core_profiles.profiles_1d[])

    @testset "simple_trust with forward_ad" begin
        dd_test = deepcopy(dd)
        act_test = deepcopy(act)
        act_test.ActorFluxMatcher.algorithm = :simple_trust

        actor = FUSE.ActorFluxMatcher(dd_test, act_test)

        # Verify the actor ran without error and produced transport results
        @test !isempty(dd_test.core_transport.model)
        cp1d_after = dd_test.core_profiles.profiles_1d[]
        @test !all(cp1d_after.electrons.temperature .== cp1d_before.electrons.temperature)
    end

    @testset "AD vs finite-diff Jacobian agreement" begin
        dd_test = deepcopy(dd)
        act_ad = deepcopy(act)
        act_ad.ActorFluxMatcher.algorithm = :simple_trust
        act_ad.ActorFluxMatcher.jacobian_method = :forward_ad
        act_ad.ActorFluxMatcher.max_iterations = 1

        dd_fd = deepcopy(dd)
        act_fd = deepcopy(act)
        act_fd.ActorFluxMatcher.algorithm = :simple_trust
        act_fd.ActorFluxMatcher.jacobian_method = :finite_diff
        act_fd.ActorFluxMatcher.max_iterations = 1

        FUSE.ActorFluxMatcher(dd_test, act_ad)
        FUSE.ActorFluxMatcher(dd_fd, act_fd)

        # After 1 iteration from same initial condition, both should produce similar profiles
        Te_ad = dd_test.core_profiles.profiles_1d[].electrons.temperature
        Te_fd = dd_fd.core_profiles.profiles_1d[].electrons.temperature
        @test isapprox(Te_ad, Te_fd; rtol=0.05)
    end

    @testset "TGLFNN forward_ad" begin
        dd_nn = deepcopy(dd)
        act_nn = deepcopy(act)
        act_nn.ActorTGLF.model = :TGLFNN
        act_nn.ActorFluxMatcher.jacobian_method = :forward_ad
        act_nn.ActorFluxMatcher.algorithm = :basic_polyalg
        act_nn.ActorFluxMatcher.max_iterations = 3

        actor = FUSE.ActorFluxMatcher(dd_nn, act_nn)

        @test !isempty(dd_nn.core_transport.model)
        cp1d_after = dd_nn.core_profiles.profiles_1d[]
        @test !all(cp1d_after.electrons.temperature .== cp1d_before.electrons.temperature)
    end

    @testset "GKNN forward_ad" begin
        dd_gk = deepcopy(dd)
        act_gk = deepcopy(act)
        act_gk.ActorTGLF.model = :GKNN
        act_gk.ActorFluxMatcher.jacobian_method = :forward_ad
        act_gk.ActorFluxMatcher.algorithm = :basic_polyalg
        act_gk.ActorFluxMatcher.max_iterations = 3

        actor = FUSE.ActorFluxMatcher(dd_gk, act_gk)

        @test !isempty(dd_gk.core_transport.model)
        cp1d_after = dd_gk.core_profiles.profiles_1d[]
        @test !all(cp1d_after.electrons.temperature .== cp1d_before.electrons.temperature)
    end
end
