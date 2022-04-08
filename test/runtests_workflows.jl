using Revise
using FUSE
using Test

@testset "init" begin

    @testset "ITER_ods" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:ITER; init_from=:ods)
        FUSE.init(dd, par)
    end

    @testset "ITER_scalars" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:ITER; init_from=:scalars)
        FUSE.init(dd, par)
    end

    @testset "CAT" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:CAT)
        FUSE.init(dd, par)
    end

    @testset "FPP_gasc" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:FPP, init_from=:gasc)
        FUSE.init(dd, par)
    end

    @testset "FPP_ods" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:FPP, init_from=:ods)
        FUSE.init(dd, par)
    end

end