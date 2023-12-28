using FUSE
using Test

@testset "warmup" begin
    for round in (1, 2)
        dd = IMAS.dd()
        FUSE.warmup(dd)
    end
end

@testset "fluxmatcher" begin
    dd, ini, act = FUSE.init(:ITER, init_from=:scalars)
    act.ActorFluxMatcher.max_iterations = 3
    act.ActorFluxMatcher.evolve_pedestal = true
    act.ActorFluxMatcher.evolve_densities = :flux_match
    FUSE.ActorFluxMatcher(dd, act)
end

@testset "init" begin
    for (testname, (args, kw)) in FUSE.test_cases
        @testset "$testname" begin
            FUSE.TimerOutputs.reset_timer!(FUSE.timer, testname)
            FUSE.TimerOutputs.@timeit FUSE.timer "$testname" begin
                println("== $(testname) ==")
                ini, act = FUSE.case_parameters(args...; kw...)
                FUSE.init(ini, act)
            end
        end
    end
end

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
