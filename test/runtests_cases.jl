using FUSE
using Test

@testset "use_cases" begin
    for (testname, (args, kw)) in FUSE.test_cases
        @testset "$testname" begin
            FUSE.TimerOutputs.reset_timer!(FUSE.timer, testname)
            FUSE.TimerOutputs.@timeit FUSE.timer "$testname" begin
                println("== $(testname) ==")

                ini, act = FUSE.case_parameters(args...; kw...)
                dd = IMAS.dd()

                # speedup the tests
                act.ActorStationaryPlasma.max_iter = 2
                # use full model for ActorThermalPlant
                act.ActorThermalPlant.model = :network

                @testset "init" begin
                    FUSE.init(dd, ini, act)
                end

                @testset "sol" begin
                    IMAS.sol(dd)
                end

                @testset "whole_facility" begin
                    FUSE.ActorWholeFacility(dd, act)
                end

                @testset "ini_dict" begin
                    FUSE.dict2ini(FUSE.ini2dict(ini))
                end

                @testset "act_dict" begin
                    FUSE.dict2act(FUSE.act2dict(act))
                end

                @testset "ini_yaml" begin
                    ini.general.dd = missing # general.dd cannot be serialized
                    ini_str = FUSE.SimulationParameters.par2ystr(ini; skip_defaults=true, show_info=false)
                    ini2 = FUSE.SimulationParameters.ystr2par(ini_str, FUSE.ParametersInits())
                end

                @testset "act_yaml" begin
                    act_str = FUSE.SimulationParameters.par2ystr(act; skip_defaults=true, show_info=false)
                    act2 = FUSE.SimulationParameters.ystr2par(act_str, FUSE.ParametersActors())
                end

                @testset "ini_json" begin
                    ini.general.dd = missing # general.dd cannot be serialized
                    ini_str = FUSE.SimulationParameters.par2jstr(ini)
                    ini2 = FUSE.SimulationParameters.jstr2par(ini_str, FUSE.ParametersInits())
                end

                @testset "act_json" begin
                    act_str = FUSE.SimulationParameters.par2jstr(act)
                    act2 = FUSE.SimulationParameters.jstr2par(act_str, FUSE.ParametersActors())
                end
            end
        end
    end
end
