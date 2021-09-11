@testset "parameters_set_get" begin
    # begin adding parameters to simple dictionary
    run_parameters = Dict()

    # scalar parameter
    run_parameters[:testScalar] = FUSE.ScalarParameter(0.5, "m/s", "scalar test")
    
    # switch parameter
    options = Dict()
    options[:bla] = FUSE.SwitchOption(:bla, "m/s", "bla")
    options[:ta] = FUSE.SwitchOption(:ta, "1/s", "ta")
    options[:user] = FUSE.ScalarParameter(1.2, "1/s", "user defined")
    run_parameters[:testSwitch] = FUSE.SwitchParameter(options, :ta, "switch test")

    # test that defining a default that is not an option throws an error 
    @test_throws Exception FUSE.SwitchParameter(options, :does_not_exist, "failing switch test")
    
    # define FuseParameters for rapid access of .value property
    local_parameters = FUSE.FuseParameters(run_parameters)
    
    @test local_parameters[:testScalar] == 0.5
    # test mutability
    local_parameters[:testScalar] = 1.0
    @test local_parameters[:testScalar] == 1.0

    @test local_parameters[:testSwitch] == :ta
    # test mutability
    local_parameters[:testSwitch] = :bla
    @test local_parameters[:testSwitch] == :bla
    # test setting something that is not a valid option
    @test_throws Exception local_parameters[:testSwitch] = :does_not_exist

    # test user defined value
    local_parameters[:testSwitch] = :user
    @test local_parameters[:testSwitch] == 1.2
    local_parameters[:testSwitch] = :user => 1.0
    @test local_parameters[:testSwitch] == 1.0
    local_parameters[:testSwitch] = :user => 2.0
    @test local_parameters[:testSwitch] == 2.0
    @test_throws Exception local_parameters[:testSwitch] = 3.0
    @test_throws Exception local_parameters[:testSwitch] = :does_not_exist => 1.0

    # test access to some parameter attribute
    @test local_parameters.parameters[:testSwitch].description == "switch test"
end

@testset "parameters_access" begin
    # test the fuse_parameters
    @test typeof(FUSE.fuse_parameters[:PHYSICS_MODELS][:bootstrapModel]) <: FUSE.CBSGi
    FUSE.fuse_parameters[:PHYSICS_MODELS][:bootstrapModel] = :user => 0.5
    @test FUSE.fuse_parameters[:PHYSICS_MODELS][:bootstrapModel] == 0.5
end