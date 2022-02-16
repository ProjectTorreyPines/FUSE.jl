module FUSE

__precompile__(true)

using IMAS
using RecipesBase
using Printf

#= ======= =#
#  INCLUDE  #
#= ======= =#
include("utils.jl")
include("physics.jl")
include("parameters.jl")

#= ====== =#
#  DDINIT  #
#= ====== =#
include("ddinit/init_equilibrium.jl")
include("ddinit/init_build.jl")
include("ddinit/init_core_profiles.jl")
include("ddinit/init_core_sources.jl")
include("ddinit/init_pf_active.jl")

#= ====== =#
#  ACTORS  #
#= ====== =#
abstract type AbstractActor end
function finalize(actor::AbstractActor)
    actor
end
include("actors/equilibrium_actor.jl")
include("actors/coils_actor.jl")
include("actors/build_actor.jl")
include("actors/current_actor.jl")
include("actors/transport_actor.jl")

#= ========= =#
#  WORKFLOWS  #
#= ========= =#
include("workflows.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export init, step, finalize
export IMAS, is_missing, @ddtime, @coords, constants, ±
export fuse_parameters, plasma_parameters, physics_models
export Parameters

end
