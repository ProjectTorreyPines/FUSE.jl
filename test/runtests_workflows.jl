using FUSE
using Test

@testset "warmup" begin
    for round in (1, 2)
        FUSE.warmup()
    end
end

@testset "init" begin
    tests = FUSE.OrderedCollections.OrderedDict()
    tests["ITER_ods"] = ([:ITER], Dict(:init_from => :ods))
    tests["ITER_scalars"] = ([:ITER], Dict(:init_from => :scalars))
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
                if occursin("ods", testname)
                    act.ActorEquilibrium.model = :Solovev
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

# @testset "QEDcurrent_actor" begin
#     # Load TRANSP data at 2.91 s
#     file_0 = joinpath(dirname(dirname(@__DIR__)), "QED", "sample", "ods_163303Z29-2910.json")
#     dd = IMAS.json2imas(file_0; verbose=false)
#     # initialize actor
#     actor = FUSE.ActorQEDcurrent(dd)
#     # evolve current
#     for k in 1:3
#         FUSE.step(actor, 0.1, 100, resume=true)
#         FUSE.finalize(actor)
#         dd.global_time = dd.equilibrium.time[end]
#     end
# end