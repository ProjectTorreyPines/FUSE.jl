using Revise
using FUSE
using Test

@testset "parameters" begin
    par = FUSE.Parameters()

    @test typeof(par.equilibrium) <: FUSE.Parameters

    @test_throws FUSE.NotsetParameterException par.equilibrium.B0

    par.equilibrium.B0 = 1.0
    @test par.equilibrium.B0 == 1.0

    @test_throws Exception par.equilibrium.B0 = "a string"

    par.equilibrium.B0 = missing
    @test_throws FUSE.NotsetParameterException par.equilibrium.B0

    @test_throws FUSE.InexistentParameterException par.equilibrium.does_not_exist = 1.0

    @test fieldnames(par.equilibrium) == fieldnames(FUSE.Parameters(:equilibrium))

    @test_throws FUSE.InexistentParameterException FUSE.Parameters(:does_not_exist)

    iter = FUSE.Parameters(:ITER; init_from=:scalars)

    @test iter.equilibrium.B0 == 5.3

    @test (iter.general.init_from = :ods) == :ods

    @test_throws FUSE.BadParameterException iter.general.init_from = :odsa
end
