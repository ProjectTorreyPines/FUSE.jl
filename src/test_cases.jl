using Test

function ini_act_tests_customizations!(ini::ParametersAllInits, act::ParametersAllActors)
    # speedup the tests
    act.ActorStationaryPlasma.max_iterations = 2
    # use full model for ActorThermalPlant if environmental variable `FUSE_WITH_EXTENSIONS` is set
    if get(ENV, "FUSE_WITH_EXTENSIONS", "false") == "true"
        act.ActorThermalPlant.model = :network
    end
    return (ini=ini, act=act)
end

function test_case(::Val{:ITER_ods}, dd::IMAS.DD)
    ini, act = case_parameters(:ITER; init_from=:ods)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:ITER_scalars}, dd::IMAS.DD)
    ini, act = case_parameters(:ITER; init_from=:scalars)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:D3D_Hmode}, dd::IMAS.DD)
    ini, act = case_parameters(:D3D, :H_mode)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:D3D_Lmode}, dd::IMAS.DD)
    ini, act = case_parameters(:D3D, :L_mode)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:D3D}, dd::IMAS.DD)
    ini, act = case_parameters(:D3D, :default)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:FPP}, dd::IMAS.DD)
    ini, act = case_parameters(:FPP)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:CAT}, dd::IMAS.DD)
    ini, act = case_parameters(:CAT)
    ini_act_tests_customizations!(ini, act)
    act.ActorStationaryPlasma.max_iterations = 1
    # There could be an engineering problem here but this shouldn't fail the tests
    act.ActorHFSsizing.error_on_performance = false
    act.ActorHFSsizing.error_on_technology = false

    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:JET_HDB5}, dd::IMAS.DD)
    ini, act = case_parameters(:HDB5; tokamak=:JET, database_case=500)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:ARC}, dd::IMAS.DD)
    ini, act = case_parameters(:ARC)
    # There could be an engineering problem here but this shouldn't fail the tests
    act.ActorHFSsizing.error_on_performance = false
    act.ActorHFSsizing.error_on_technology = false

    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:SPARC}, dd::IMAS.DD)
    ini, act = case_parameters(:SPARC; init_from=:ods)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:KDEMO}, dd::IMAS.DD)
    ini, act = case_parameters(:KDEMO)
    ini_act_tests_customizations!(ini, act)
    act.ActorFluxMatcher.max_iterations = 5
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:KDEMO_compact}, dd::IMAS.DD)
    ini, act = case_parameters(:KDEMO_compact)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:DTT}, dd::IMAS.DD)
    ini, act = case_parameters(:DTT)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:EXCITE}, dd::IMAS.DD)
    ini, act = case_parameters(:EXCITE)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:MANTA}, dd::IMAS.DD)
    ini, act = case_parameters(:MANTA)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:UNIT}, dd::IMAS.DD)
    ini, act = case_parameters(:UNIT)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:ITER_time}, dd::IMAS.DD; verbose::Bool=false, Δt::Float64=60.0, Nt::Int=12)
    # ========
    # hardware setup from ODS

    ini, act = FUSE.case_parameters(:ITER; init_from=:ods);
    FUSE.init(dd, ini, act);

    # ========
    # pulse_schedule fom scalars

    ini, _ = FUSE.case_parameters(:ITER; init_from=:scalars, time_dependent=true);
    ini.time.simulation_start = 50.0
    FUSE.init(dd, ini, act; initialize_hardware=false);

    # ========
    # Our simulation should start in a self-consistent state. For this, we call the `ActorStationaryPlasma`

    act.ActorStationaryPlasma.convergence_error = 2E-2
    act.ActorStationaryPlasma.max_iterations = 10

    act.ActorFluxMatcher.verbose = verbose
    act.ActorFluxMatcher.relax = 1.0
    act.ActorFluxMatcher.algorithm = :simple_dfsane
    act.ActorSawteethSource.flat_factor = 1.0

    FUSE.ActorStationaryPlasma(dd, act; verbose)

    # ========
    # Now we're ready to actually run the time-dependent simulation

    # defaults evolve 60s at the same 5s time step as the full 300s simulation,
    # to keep this test short: test_case(:ITER_time; Δt=300.0, Nt=60) for the full run
    act.ActorDynamicPlasma.Nt = Nt
    act.ActorDynamicPlasma.Δt = Δt

    act.ActorDynamicPlasma.evolve_current = true
    act.ActorDynamicPlasma.evolve_equilibrium = true
    act.ActorDynamicPlasma.evolve_transport = true
    act.ActorDynamicPlasma.evolve_sources = true
    act.ActorDynamicPlasma.evolve_pf_active = false
    act.ActorDynamicPlasma.evolve_pedestal = true
    act.ActorDynamicPlasma.evolve_sawteeth = true

    act.ActorDynamicPlasma.ip_controller = true
    act.ActorDynamicPlasma.time_derivatives_sources = true
    FUSE.ActorDynamicPlasma(dd, act; verbose);

    return (dd=dd, ini=ini, act=act)
end

# ================ #

function test_case(case::Symbol)
    dd = IMAS.dd()
    return test_case(Val(case), dd::IMAS.DD)
end

function test_case(case::Symbol, dd::IMAS.DD)
    dd, ini, act = test_case(Val(case), dd::IMAS.DD)
    @assert typeof(dd) <: IMAS.dd
    @assert typeof(ini) <: ParametersAllInits
    @assert typeof(act) <: ParametersAllActors
    return (dd=dd, ini=ini, act=act)
end

function available_test_cases()
    return sort!([m.sig.parameters[2].parameters[1] for m in methods(test_case) if m.sig.parameters[2] <: Val])
end

# ================ #

function test_ini_act_save_load(dd::IMAS.DD, ini::ParametersAllInits, act::ParametersAllActors)
    Test.@testset "init" begin
        init(dd, ini, act)
    end

    Test.@testset "sol" begin
        IMAS.sol(dd)
    end

    Test.@testset "whole_facility" begin
        ActorWholeFacility(dd, act)
    end

    Test.@testset "ini_dict" begin
        dict2ini(ini2dict(ini))
    end

    Test.@testset "act_dict" begin
        dict2act(act2dict(act))
    end

    Test.@testset "ini_yaml" begin
        ini.general.dd = missing
        ini_str = SimulationParameters.par2ystr(ini; skip_defaults=true, show_info=false)
        ini2 = SimulationParameters.ystr2par(ini_str, ParametersInits())
    end

    Test.@testset "act_yaml" begin
        act_str = SimulationParameters.par2ystr(act; skip_defaults=true, show_info=false)
        act2 = SimulationParameters.ystr2par(act_str, ParametersActors())
    end

    Test.@testset "ini_json" begin
        ini.general.dd = missing
        ini_str = SimulationParameters.par2jstr(ini)
        ini2 = SimulationParameters.jstr2par(ini_str, ParametersInits())
    end

    Test.@testset "act_json" begin
        act_str = SimulationParameters.par2jstr(act)
        act2 = SimulationParameters.jstr2par(act_str, ParametersActors())
    end

    Test.@testset "ini_hdf5" begin
        tmpdir = mktempdir()
        ini.general.dd = missing
        SimulationParameters.par2hdf(ini, joinpath(tmpdir, "ini.h5"))
        ini2 = SimulationParameters.hdf2par(joinpath(tmpdir, "ini.h5"), ParametersInits())
        rm(tmpdir; force=true, recursive=true)
    end

    Test.@testset "act_hdf5" begin
        tmpdir = mktempdir()
        SimulationParameters.par2hdf(act, joinpath(tmpdir, "act.h5"))
        act2 = SimulationParameters.hdf2par(joinpath(tmpdir, "act.h5"), ParametersActors())
        rm(tmpdir; force=true, recursive=true)
    end

    return (ini=ini, act=act)
end
