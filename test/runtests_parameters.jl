using Revise
using FUSE
using Test
using InteractiveUtils: subtypes

@testset "ParametersInit" begin
    par = FUSE.ParametersInit()

    @test typeof(par.equilibrium) <: FUSE.ParametersInit

    @test_throws FUSE.NotsetParameterException par.equilibrium.B0

    par.equilibrium.B0 = 1.0
    @test par.equilibrium.B0 == 1.0

    @test par.equilibrium[:B0]._name == :B0
    @test par[:equilibrium]._name == :equilibrium

    ini = FUSE.ParametersInit()
    ini1 = FUSE.ParametersInit()
    ini.tf = ini1.tf
    @test ini.tf._parent.value !== ini1.tf._parent.value

    @test_throws Exception par.equilibrium.B0 = "a string"

    par.equilibrium.B0 = missing
    @test_throws FUSE.NotsetParameterException par.equilibrium.B0

    @test_throws FUSE.InexistentParameterException par.equilibrium.does_not_exist = 1.0

    @test keys(par.equilibrium) == keys(FUSE.ParametersInit(:equilibrium))

    @test_throws FUSE.InexistentParameterException FUSE.ParametersInit(:does_not_exist)

    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)

    @test typeof(ini) <: FUSE.ParametersInit

    @test typeof(act) <: FUSE.ParametersActor

    @test ini.equilibrium.B0 == -5.3

    @test (ini.general.init_from = :ods) == :ods

    @test_throws FUSE.BadParameterException ini.general.init_from = :odsa

    @test_throws FUSE.InexistentParameterException FUSE.ParametersInit(:inexistent_group)

    @test_throws UndefKeywordError FUSE.case_parameters(:ITER)

    for par in subtypes(FUSE.AbstractActor)
        par = Symbol(replace(string(par), "FUSE." => ""))
        @test typeof(FUSE.ParametersActor(par)) <: FUSE.ParametersActor
    end

end
