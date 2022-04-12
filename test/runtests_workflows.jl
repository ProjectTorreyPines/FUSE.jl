using Revise
using FUSE
using Test

@testset "init" begin

    @testset "ITER_ods" begin
        dd, ini, act = FUSE.init(:ITER; init_from=:ods)
    end

    @testset "ITER_scalars" begin
        dd, ini, act = FUSE.init(:ITER; init_from=:scalars)
    end

    @testset "D3D" begin
        dd, ini, act = FUSE.init(:D3D)
    end

    @testset "FPP_ods" begin
        dd, ini, act = FUSE.init(:FPP; init_from=:ods)
    end

    @testset "FPP_gasc" begin
        dd, ini, act = FUSE.init(:FPP; init_from=:gasc)
    end

    @testset "CAT" begin
        dd, ini, act = FUSE.init(:CAT)
    end

    @testset "HDB5" begin
        dd, ini, act = FUSE.init(:HDB5; tokamak=:JET, case=500)
    end

end


# @testset "QEDcurrent_actor" begin
#     # Load TRANSP data at 2.91 s
#     file_0 = joinpath(dirname(dirname(dirname(abspath(@__FILE__)))), "QED", "sample", "ods_163303Z29-2910.json")
#     dd = IMAS.json2imas(file_0; verbose=false)
#     # initialize actor
#     actor = FUSE.QEDcurrentActor(dd)
#     # evolve current
#     for k in 1:3
#         FUSE.step(actor, 0.1, 100, resume=true)
#         FUSE.finalize(actor)
#         dd.global_time = dd.equilibrium.time[end]
#     end
# end