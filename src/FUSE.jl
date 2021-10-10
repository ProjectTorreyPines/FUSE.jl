__precompile__()

module FUSE

using IMAS

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
export IMAS
export fuse_parameters, plasma_parameters, physics_models

end # module
