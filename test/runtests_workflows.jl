using FUSE
using Test

@testset "init" begin

    @testset "ITER_ods" begin
        dd, ini, act = FUSE.init(:ITER; init_from=:ods)
    end

    @testset "ITER_scalars" begin
        dd, ini, act = FUSE.init(:ITER; init_from=:scalars)
    end

    @testset "D3D" begin
        dd, ini, act = FUSE.init(:D3D)
    end

    @testset "FPP_v1_demount_scalars" begin
        dd, ini, act = FUSE.init(:FPP; version=:v1_demount, init_from=:scalars)
    end

    @testset "FPP_v1_demount_ods" begin
        dd, ini, act = FUSE.init(:FPP; version=:v1_demount, init_from=:ods)
    end

    @testset "FPP_v1_ods" begin
        dd, ini, act = FUSE.init(:FPP; version=:v1, init_from=:ods)
    end

    @testset "FPP_v1_scalars" begin
        dd, ini, act = FUSE.init(:FPP; version=:v1, init_from=:scalars)
    end

    @testset "CAT" begin
        dd, ini, act = FUSE.init(:CAT)
    end

    @testset "HDB5" begin
        dd, ini, act = FUSE.init(:HDB5; tokamak=:JET, case=500)
    end

    @testset "ARC" begin
        dd, ini, act = FUSE.init(:ARC)
    end

    @testset "SPARC" begin
        dd, ini, act = FUSE.init(:SPARC)
    end
end

@testset "warmup" begin
    FUSE.warmup()
end

@testset "optimization" begin
    ini = FUSE.ParametersInits()

    ini.core_profiles.zeff = 2.0 ↔ [1.2, 2.5]
    @test ini.core_profiles.zeff == 2.0
    @test typeof(ini.core_profiles.zeff) <: Float64

    ini.equilibrium.ngrid = 200.0 ↔ [129, 257]
    @test ini.equilibrium.ngrid == 200
    @test typeof(ini.equilibrium.ngrid) <: Int

    ini.build.symmetric = 1.0 ↔ [0, 1]
    @test ini.build.symmetric == true
    @test typeof(ini.build.symmetric) <: Bool
end

# @testset "QEDcurrent_actor" begin
#     # Load TRANSP data at 2.91 s
#     file_0 = joinpath(dirname(dirname(@__DIR__)), "QED", "sample", "ods_163303Z29-2910.json")
#     dd = IMAS.json2imas(file_0; verbose=false)
#     # initialize actor
#     actor = FUSE.ActorQEDcurrent(dd)
#     # evolve current
#     for k in 1:3
#         FUSE.step(actor, 0.1, 100, resume=true)
#         FUSE.finalize(actor)
#         dd.global_time = dd.equilibrium.time[end]
#     end
# end