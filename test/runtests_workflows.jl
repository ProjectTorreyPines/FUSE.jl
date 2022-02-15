using Revise
using FUSE
using Test

@testset "workflows" begin

    @testset "ITER" begin
        par = FUSE.Parameters(:ITER)
        par.general.init_from = :ods
        FUSE.simple_workflow(par)

        par = FUSE.Parameters(:ITER)
        par.general.init_from = :scalars
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