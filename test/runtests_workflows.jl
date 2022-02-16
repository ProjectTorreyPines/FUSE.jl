using Revise
using FUSE
using Test

@testset "workflows" begin

    @testset "ITER_ods" begin
        par = FUSE.Parameters(:ITER; init_from = :ods)
        FUSE.simple_workflow(par)
    end

    @testset "ITER_scalars" begin
        par = FUSE.Parameters(:ITER; init_from = :scalars)
        FUSE.simple_workflow(par)
    end

    @testset "CAT" begin
        par = FUSE.Parameters(:CAT)
        FUSE.simple_workflow(par)
    end

    @testset "FPP" begin
        par = FUSE.Parameters(:FPP)
        FUSE.simple_workflow(par)
    end

end