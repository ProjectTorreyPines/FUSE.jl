__precompile__()

module FUSE

import IMAS

#= ========= =#
#= BOOTSTRAP =#
#= ========= =#

include("bootstrap.jl")

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
export IMAS
export fuse_parameters, plasma_parameters, physics_models

end # module
