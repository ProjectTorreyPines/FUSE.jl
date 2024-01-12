using FUSE
using Test

@testset "optimization" begin
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
