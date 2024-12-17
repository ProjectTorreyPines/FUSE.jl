using FUSE
using Test

@testset "warmup" begin
    println("== warmup ==")
    for round in (1, 2)
        dd = IMAS.dd()
        @test begin
            FUSE.warmup(dd)
            true
        end
    end
end