using Revise
using FUSE
using Test



@testset "parameters_set_get" begin
    PAR = AllParameters()

    @test typeof(PAR.equilibrium) <: FUSE.Parameters

    @test_throws FUSE.NotsetParameterException PAR.equilibrium.B0

    PAR.equilibrium.B0 = 1.0
    @test PAR.equilibrium.B0 == 1.0

    @test_throws Exception PAR.equilibrium.B0 = "a string"

    PAR.equilibrium.B0 = missing
    @test_throws FUSE.NotsetParameterException PAR.equilibrium.B0

    @test_throws FUSE.InexistentParameterException PAR.equilibrium.does_not_exist = 1.0

end
