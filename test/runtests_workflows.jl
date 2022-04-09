using Revise
using FUSE
using Test

@testset "init" begin

    @testset "ITER_ods" begin
        dd = IMAS.dd()
        par = FUSE.InitParameters(:ITER; init_from=:ods)
        FUSE.init(dd, par)
    end

    @testset "ITER_scalars" begin
        dd = IMAS.dd()
        par = FUSE.InitParameters(:ITER; init_from=:scalars)
        FUSE.init(dd, par)
    end

    @testset "CAT" begin
        dd = IMAS.dd()
        par = FUSE.InitParameters(:CAT)
        FUSE.init(dd, par)
    end

    @testset "FPP_gasc" begin
        dd = IMAS.dd()
        par = FUSE.InitParameters(:FPP, init_from=:gasc)
        FUSE.init(dd, par)
    end

    @testset "FPP_ods" begin
        dd = IMAS.dd()
        par = FUSE.InitParameters(:FPP, init_from=:ods)
        FUSE.init(dd, par)
    end

    @testset "HDB5" begin
        dd = IMAS.dd()
        par = FUSE.InitParameters(:HDB5; tokamak=:JET, case=500)
        FUSE.init(dd, par)
    end

end


@testset "QEDcurrent_actor" begin
    # Load TRANSP data at 2.91 s
    file_0 = joinpath(dirname(dirname(dirname(abspath(@__FILE__)))), "QED", "sample", "ods_163303Z29-2910.json")
    dd = IMAS.json2imas(file_0; verbose=false)
    # initialize actor
    actor = FUSE.QEDcurrentActor(dd)
    # evolve current
    for k in 1:3
        FUSE.step(actor, 0.1, 100, resume=true)
        FUSE.finalize(actor)
        dd.global_time = dd.equilibrium.time[end]
    end
end