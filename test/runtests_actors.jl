using FUSE
using Test

@testset "ActorGPEC" begin
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd = FUSE.init(ini, act)

    @test length(dd.equilibrium.time_slice) > 0

    act.ActorGPEC.nn_low   = 1
    act.ActorGPEC.nn_high  = 1
    act.ActorGPEC.vac_flag = true
    act.ActorGPEC.verbose  = false
    # Keep psihigh below GPEC's internal edge-scan threshold (psiedge=0.99,
    # not exposed via the actor) so the run stays on the well-tested core path.
    act.ActorGPEC.psihigh  = 0.95

    FUSE.ActorGPEC(dd, act)

    @test length(dd.mhd_linear.time_slice) > 0
    @test length(dd.mhd_linear.time_slice[1].toroidal_mode) == 1
    @test real(dd.mhd_linear.time_slice[1].toroidal_mode[1].energy_perturbed) > 0
end

@testset "fluxmatcher" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    act.ActorFluxMatcher.max_iterations = 2
    act.ActorFluxMatcher.evolve_pedestal = true

    act.ActorFluxMatcher.algorithm = :simple
    act.ActorTGLF.model = :TJLF
    FUSE.ActorFluxMatcher(dd, act)
    act.ActorFluxMatcher.algorithm = :anderson
    act.ActorTGLF.model = :TGLFNN
    FUSE.ActorFluxMatcher(dd, act)
end
