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
    act.ActorFluxMatcher.evolve_densities = FUSE.setup_density_evolution_electron_flux_match_rest_ne_scale(dd)
    FUSE.ActorFluxMatcher(dd, act)
end

@testset "dd floats" begin
    ini, act = FUSE.init(:ITER, init_from=:ods)
    dd = IMAS.dd{Real}()
    FUSE.init(dd, ini, act)
end

@testset "init" begin
    tests = FUSE.OrderedCollections.OrderedDict()
    tests["D3D"] = ([:D3D], Dict())
    tests["FPP_v1_demount_scalars"] = ([:FPP], Dict(:version => :v1_demount, :init_from => :scalars))
    tests["FPP_v1_demount_ods"] = ([:FPP], Dict(:version => :v1_demount, :init_from => :ods))
    tests["FPP_v1_ods"] = ([:FPP], Dict(:version => :v1, :init_from => :ods))
    tests["FPP_v1_scalars"] = ([:FPP], Dict(:version => :v1, :init_from => :scalars))
    tests["CAT"] = ([:CAT], Dict())
    tests["HDB5"] = ([:HDB5], Dict(:tokamak => :JET, :case => 500))
    tests["ARC"] = ([:ARC], Dict())
    tests["SPARC"] = ([:SPARC], Dict())

    for (testname, (args, kw)) in tests
        @testset "$testname" begin
            FUSE.TimerOutputs.reset_timer!(FUSE.timer, testname)
            FUSE.TimerOutputs.@timeit FUSE.timer "$testname" begin
                println("== $(testname) ==")
                ini, act = FUSE.case_parameters(args...; kw...)
                if testname == "SPARC"
                    act.ActorEquilibrium.model = :Solovev # The CI tool doesn't work with SPARC & TEQUILA for unexplainable reasons
                end
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
