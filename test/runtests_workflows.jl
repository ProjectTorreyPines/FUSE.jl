using Revise
using FUSE
using Test

@testset "workflows" begin

    @testset "ITER_ods" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:ITER; init_from = :ods)
        FUSE.init_workflow(dd, par)
    end

    @testset "ITER_scalars" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:ITER; init_from = :scalars)
        FUSE.init_workflow(dd, par)
    end

    @testset "CAT" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:CAT)
        FUSE.init_workflow(dd, par)
    end

    @testset "FPP" begin
        dd = IMAS.dd()
        par = FUSE.Parameters(:FPP)
        FUSE.init_workflow(dd, par)
    end

end