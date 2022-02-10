__precompile__()

module FUSE

using IMAS
using Printf

#= ======= =#
#= INCLUDE =#
#= ======= =#
include("utils.jl")

include("physics.jl")

include("actors.jl")

include("parameters.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export init, step, finalize
export IMAS, is_missing, @ddtime, @coords, constants, ±
export fuse_parameters, plasma_parameters, physics_models

end
