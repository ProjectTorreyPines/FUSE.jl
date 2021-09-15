__precompile__()

module FUSE

using IMAS

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
# export f2p, f2u, p2i, i2p, u2f, f2f, f2fs, u2fs
# export coordinates
# export IDS, IDSvector, IDSvectorElement
# export dd
# export top, parent, children, expressions
export IMAS, set_field_time_array
export fuse_parameters, plasma_parameters, physics_models

end # module
