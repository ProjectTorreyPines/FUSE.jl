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
                act.ActorStationaryPlasma.max_iter = 1

                @testset "init" begin
                    FUSE.init(dd, ini, act)
                end

                @testset "sol" begin
                    IMAS.sol(dd)
                end

                @testset "whole_facility" begin
                    FUSE.ActorWholeFacility(dd, act)
                end

                @testset "ini_yaml" begin
                    ini.general.dd = missing
                    ini_str = FUSE.SimulationParameters.par2ystr(ini; skip_defaults=true, show_info=false)
                    ini2 = try
                        FUSE.SimulationParameters.ystr2par(ini_str, FUSE.ParametersInits())
                    catch e
                        rethrow(e)
                    end
                end

                @testset "act_yaml" begin
                    act_str = FUSE.SimulationParameters.par2ystr(act; skip_defaults=true, show_info=false)
                    act2 = try
                        FUSE.SimulationParameters.ystr2par(act_str, FUSE.ParametersActors())
                    catch e
                        rethrow(e)
                    end
                end

                @testset "ini_json" begin
                    ini.general.dd = missing
                    ini_str = FUSE.SimulationParameters.par2jstr(ini)
                    ini2 = try
                        FUSE.SimulationParameters.jstr2par(ini_str, FUSE.ParametersInits())
                    catch e
                        rethrow(e)
                    end
                end

                @testset "act_json" begin
                    act_str = FUSE.SimulationParameters.par2jstr(act)
                    act2 = try
                        FUSE.SimulationParameters.jstr2par(act_str, FUSE.ParametersActors())
                    catch e
                        rethrow(e)
                    end
                end
            end
        end
    end
end
