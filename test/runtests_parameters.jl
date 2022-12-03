using Revise
using FUSE
using Test

@testset "ParametersInit" begin
    par = FUSE.ParametersAllInits()

    @test typeof(par.equilibrium) <: FUSE.AbstractParameters

    @test_throws FUSE.NotsetParameterException par.equilibrium.B0

    par.equilibrium.B0 = 1.0
    @test par.equilibrium.B0 == 1.0

    @test par.equilibrium[:B0]._name == :B0
    @test getfield(par[:equilibrium],:_name) == :equilibrium

    ini = FUSE.ParametersAllInits()
    ini1 = FUSE.ParametersAllInits()
    ini.tf = ini1.tf
    @test getfield(ini.tf, :_parent).value !== getfield(ini1.tf, :_parent).value

    @test_throws Exception par.equilibrium.B0 = "a string"

    par.equilibrium.B0 = missing
    @test_throws FUSE.NotsetParameterException par.equilibrium.B0

    @test_throws FUSE.InexistentParameterException par.equilibrium.does_not_exist = 1.0

    @test keys(par.equilibrium) == keys(FUSE.ParametersInit(:equilibrium))

    @test_throws FUSE.InexistentParameterException FUSE.ParametersInit(:does_not_exist)

    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)

    @test typeof(ini) <: FUSE.AbstractParameters

    @test typeof(act) <: FUSE.AbstractParameters

    @test ini.equilibrium.B0 == -5.3

    @test (ini.general.init_from = :ods) == :ods

    @test_throws FUSE.BadParameterException ini.general.init_from = :odsa

    @test_throws FUSE.InexistentParameterException FUSE.ParametersInit(:inexistent_group)

    @test_throws UndefKeywordError FUSE.case_parameters(:ITER)

    # fail tests if not all actors have parameters associated with them
    for par in FUSE.concretetypes(FUSE.AbstractActor)
        par = Symbol(replace(string(par), "FUSE." => ""))
        @test typeof(FUSE.ParametersActor(par)) <: FUSE.ParametersActor
    end

    # save load
    FUSE.dict2par!(FUSE.par2dict(ini), FUSE.ParametersAllInits())
    FUSE.dict2par!(FUSE.par2dict(act), FUSE.ParametersAllActors())

end
