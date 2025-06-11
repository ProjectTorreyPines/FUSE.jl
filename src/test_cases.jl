using Test

function test_case(::Val{:ITER_ods}, dd::IMAS.dd)
    ini, act = case_parameters(:ITER; init_from=:ods)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:ITER_scalars}, dd::IMAS.dd)
    ini, act = case_parameters(:ITER; init_from=:scalars)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:D3D_Hmode}, dd::IMAS.dd)
    ini, act = case_parameters(:D3D, :H_mode)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:D3D_Lmode}, dd::IMAS.dd)
    ini, act = case_parameters(:D3D, :L_mode)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:D3D}, dd::IMAS.dd)
    ini, act = case_parameters(:D3D, :default)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:FPP}, dd::IMAS.dd)
    ini, act = case_parameters(:FPP)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:CAT}, dd::IMAS.dd)
    ini, act = case_parameters(:CAT)
    ini_act_tests_customizations!(ini, act)
    act.ActorStationaryPlasma.max_iterations = 1
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:JET_HDB5}, dd::IMAS.dd)
    ini, act = case_parameters(:HDB5; tokamak=:JET, database_case=500)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:ARC}, dd::IMAS.dd)
    ini, act = case_parameters(:ARC)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:SPARC}, dd::IMAS.dd)
    ini, act = case_parameters(:SPARC; init_from=:ods)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:KDEMO}, dd::IMAS.dd)
    ini, act = case_parameters(:KDEMO)
    ini_act_tests_customizations!(ini, act)
    act.ActorFluxMatcher.max_iterations = 5
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:KDEMO_compact}, dd::IMAS.dd)
    ini, act = case_parameters(:KDEMO_compact)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:DTT}, dd::IMAS.dd)
    ini, act = case_parameters(:DTT)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:EXCITE}, dd::IMAS.dd)
    ini, act = case_parameters(:EXCITE)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:MANTA}, dd::IMAS.dd)
    ini, act = case_parameters(:MANTA)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

function test_case(::Val{:UNIT}, dd::IMAS.dd)
    ini, act = case_parameters(:UNIT)
    ini_act_tests_customizations!(ini, act)
    test_ini_act_save_load(dd, ini, act)
    return (dd=dd, ini=ini, act=act)
end

# ================ #

function test_case(case::Symbol, dd::IMAS.dd)
    return test_case(Val(case), dd::IMAS.dd)
end

function available_test_cases()
    return sort!([m.sig.parameters[2].parameters[1] for m in methods(test_case) if m.sig.parameters[2] <: Val])
end

# ================ #

function ini_act_tests_customizations!(ini::ParametersAllInits, act::ParametersAllActors)
    # speedup the tests
    act.ActorStationaryPlasma.max_iterations = 2
    # use full model for ActorThermalPlant if environmental variable `FUSE_WITH_EXTENSIONS` is set
    if get(ENV, "FUSE_WITH_EXTENSIONS", "false") == "true"
        act.ActorThermalPlant.model = :network
    end
    return (ini=ini, act=act)
end

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
