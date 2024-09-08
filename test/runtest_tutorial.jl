using FUSE
using Test

@testset "tutorial" begin
    include(joinpath(@__DIR__, "..", "docs", "src", "tutorial.jl"))
end

