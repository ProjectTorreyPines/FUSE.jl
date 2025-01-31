using FUSE
using Test

@testset "warmup" begin
    for round in ("before_compile", "after_compile")
        println("== warmup_$(round) ==")
        FUSE.TimerOutputs.@timeit FUSE.timer "warmup_$(round)" begin
            dd = IMAS.dd()
            @test begin
                FUSE.warmup(dd)
                true # this is to make the test pass
            end
        end
    end
end
