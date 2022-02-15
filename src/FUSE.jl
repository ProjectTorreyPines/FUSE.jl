module FUSE

__precompile__(true)

using IMAS
using Printf

#= ======= =#
#= INCLUDE =#
#= ======= =#
include("utils.jl")

include("physics.jl")

include("parameters.jl")

include("actors.jl")

include("workflows.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export init, step, finalize
export IMAS, is_missing, @ddtime, @coords, constants, ±
export fuse_parameters, plasma_parameters, physics_models
export Parameters

end
