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

    FUSE.ActorGPEC(dd, act)

    @test length(dd.mhd_linear.time_slice) > 0
    @test length(dd.mhd_linear.time_slice[1].toroidal_mode) == 1
    @test real(dd.mhd_linear.time_slice[1].toroidal_mode[1].energy_perturbed) > 0
end
