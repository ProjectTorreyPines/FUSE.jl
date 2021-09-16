__precompile__()

module FUSE

using IMAS
using Printf
using ForwardDiff

include("utils.jl")

#= ====== =#
#= ACTORS =#
#= ====== =#

include("actors.jl")

#= ========== =#
#= PARAMETERS =#
#= ========== =#

include("parameters.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export init, step, finalize
export SolovevEquilibriumActor
export IMAS
export fuse_parameters, plasma_parameters, physics_models

end # module
