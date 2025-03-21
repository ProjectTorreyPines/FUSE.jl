using FUSE
using Test

@testset "GC" begin
    act = FUSE.ParametersActors()
    @test parent(act.ActorNeoclassical) === act
    GC.gc()
    @test parent(act.ActorNeoclassical) === act
end

@testset "optimization_parameters" begin
    ini = FUSE.ParametersInits()

    ini.core_profiles.zeff = 2.0 ↔ [1.2, 2.5]
    @test ini.core_profiles.zeff == 2.0
    @test typeof(ini.core_profiles.zeff) <: Float64

    ini.equilibrium.ngrid = 200 ↔ [129, 257]
    @test ini.equilibrium.ngrid == 200
    @test typeof(ini.equilibrium.ngrid) <: Int

    ini.build.symmetric = true ↔ (false, true)
    @test ini.build.symmetric == true
    @test typeof(ini.build.symmetric) <: Bool
end

@testset "Checkpoint" begin
    a_old = a = [1]
    b_old = b = [2]
    dd_old = dd = IMAS.dd()
    ini_old = ini = FUSE.ParametersInits()

    # checkpoint empty
    empty!(FUSE.checkpoint)

    # checkin
    @checkin :test1 a b dd

    a[1] = 100

    @test a[1] == 100
    @test a === a_old
    @test b[1] == 2
    @test b === b_old
    @test dd === dd_old

    # add to the same checking
    @checkin :test1 ini

    # partial checkout
    @checkout :test1 a b

    @test a[1] == 1
    @test a !== a_old
    @test b[1] == 2
    @test b !== b_old

    # partial checkout
    @checkout :test1 dd ini

    @test dd !== dd_old
    @test ini !== ini_old

    # checkpoint empty
    empty!(FUSE.checkpoint)
    @assert isempty(FUSE.checkpoint)
end