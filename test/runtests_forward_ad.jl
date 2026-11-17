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
        # the ITER case selects sat1_em_iter, which GKNN does not support
        act_gk.ActorTGLF.sat_rule = :sat3
        act_gk.ActorTGLF.tglfnn_model = "sat3_em_d3d_azf-1"
        act_gk.ActorFluxMatcher.jacobian_method = :forward_ad
        act_gk.ActorFluxMatcher.algorithm = :basic_polyalg
        act_gk.ActorFluxMatcher.max_iterations = 3

        actor = FUSE.ActorFluxMatcher(dd_gk, act_gk)

        @test !isempty(dd_gk.core_transport.model)
        cp1d_after = dd_gk.core_profiles.profiles_1d[]
        @test !all(cp1d_after.electrons.temperature .== cp1d_before.electrons.temperature)
    end

    @testset "forward-mode residual equals the solver's residual" begin
        # the forward-mode evaluation builds its own dd: it must see the same sources (targets, fast-ion
        # densities) as the primal residual. Without the external sources the two differed by 11 % here
        dd_test = deepcopy(dd)
        act_test = deepcopy(act)
        act_test.ActorTGLF.model = :TGLFNN
        act_test.ActorFluxMatcher.algorithm = :none
        act_test.ActorFluxMatcher.consistent_fast_pressure = true
        actor = FUSE.ActorFluxMatcher(dd_test, act_test)
        cp1d = dd_test.core_profiles.profiles_1d[]
        initial_cp1d = FUSE.cp1d_copy_primary_quantities(cp1d)
        z = FUSE.scale_z_profiles(FUSE.pack_z_profiles(cp1d, actor.par)[1])
        primal = FUSE.flux_match_errors(actor, z, initial_cp1d).errors
        F = similar(primal)
        FUSE.ad_flux_match_errors!(F, z, actor, initial_cp1d)
        @test F ≈ primal rtol = 1e-2
    end
end

@testset "forward_ad Jacobian does not replace the solver's residual" begin
    # NonlinearSolve's TrustRegion with jacobian_method=:forward_ad. The residual of the forward-mode evaluation
    # (ad_flux_match_errors!) used to leave out the external sources, so when the solver took its residual from the
    # Jacobian evaluation it stepped from a wrong residual: on this sample it stalled at the starting error (1.2)
    dd = IMAS.json2imas(joinpath(pkgdir(FUSE), "sample", "D3D_standard_Lmode.json"); show_warnings=false)
    act = FUSE.ParametersActors()
    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"
    act.ActorFluxMatcher.algorithm = :trust
    act.ActorFluxMatcher.jacobian_method = :forward_ad
    actor = FUSE.ActorFluxMatcher(dd, act)
    @test actor.error < act.ActorFluxMatcher.xtol
end

@testset "forward_ad Jacobian with a finite time step" begin
    # With a finite Δt the residual has the source of the change the flux matcher makes to the profiles over Δt,
    # whose derivative grows as 1/Δt. Without it in the forward-mode Jacobian the solvers did not converge for short
    # time steps (on this case, at Δt = 5 ms, TrustRegion ended at an error of 0.5 after 52 evaluations)
    Δt = 0.005
    dd = IMAS.json2imas(joinpath(pkgdir(FUSE), "sample", "D3D_standard_Lmode.json"); show_warnings=false)
    t0 = dd.core_profiles.time[end]
    dd.global_time = t0 + Δt
    for ids in (dd.equilibrium, dd.core_sources, dd.core_profiles)
        IMAS.new_timeslice!(ids, t0 + Δt)
    end
    cp1d = dd.core_profiles.profiles_1d[]
    cp1d.electrons.temperature = cp1d.electrons.temperature .* 1.05
    cp1d.electrons.density_thermal = cp1d.electrons.density_thermal .* 1.03
    for ion in cp1d.ion
        if !ismissing(ion, :density_thermal)
            ion.density_thermal = ion.density_thermal .* 1.03
        end
    end
    act = FUSE.ParametersActors()
    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"
    act.ActorFluxMatcher.algorithm = :trust
    act.ActorFluxMatcher.jacobian_method = :forward_ad
    act.ActorFluxMatcher.Δt = Δt
    act.ActorFluxMatcher.relax = 1.0
    actor = FUSE.ActorFluxMatcher(dd, act)
    @test actor.error < act.ActorFluxMatcher.xtol
    @test length(actor.err_history) < 20
end
