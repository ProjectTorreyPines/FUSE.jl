using FUSE
using Test

include("runtests_workflows.jl")

include("runtests_basics.jl")

println(FUSE.to)
