using FUSE
using Test

@testset "act_migrations" begin
    act = FUSE.ParametersActors()
    schema_keys = Set(string(k) for k in keys(act))
    for old_name in keys(FUSE.ACT_ACTOR_MIGRATIONS)
        @test old_name ∉ schema_keys
    end
end

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

@testset "evolve_densities :zeff" begin
    dd = IMAS.json2imas(joinpath(pkgdir(FUSE), "sample", "D3D_standard_Lmode.json"); show_warnings=false)
    cp1d = dd.core_profiles.profiles_1d[]

    # fast electrons are not among the ions that get :zeff: they are registered as :fixed
    cp1d.electrons.density_fast = zero(cp1d.electrons.density_thermal)
    evolve_densities = FUSE.setup_density_evolution_electron_flux_match_fixed_zeff(cp1d)
    @test !haskey(evolve_densities, :electrons_fast)
    @test FUSE.check_evolve_densities(cp1d, evolve_densities) === nothing
    @test evolve_densities[:electrons_fast] == :fixed

    # a fast ion can be :fixed too
    evolve_densities[:D_fast] = :fixed
    @test FUSE.check_evolve_densities(cp1d, evolve_densities) === nothing

    # every thermal ion must be :zeff
    evolve_densities[:C] = :fixed
    @test_throws AssertionError FUSE.check_evolve_densities(cp1d, evolve_densities)
end

@testset "equalize_symmetric_areas!" begin
    # coils 2 and 3 mirror each other; coil 1 sits next to coil 2, within a coil size of the mirror image of coil 3
    # coil 4 has no mirror; coils 5 and 6 mirror each other
    r = [5.1, 5.0, 5.0, 8.0, 1.0, 1.0]
    z = [1.0, 1.0, -1.0, 0.5, 2.0, -2.0]
    areas = [1.0, 2.0, 3.0, 0.5, 4.0, 6.0]
    @test FUSE.equalize_symmetric_areas!(areas, r, z) == [1.0, 3.0, 3.0, 0.5, 6.0, 6.0]
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