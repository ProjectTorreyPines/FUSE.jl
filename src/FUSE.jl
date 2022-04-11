module FUSE

__precompile__(true)

using IMAS
using Plots
using Printf

#= ======= =#
#  UTILS  #
#= ======= =#
include("utils.jl")

#= ========== =#
#  PARAMETERS  #
#= ========== =#
include("parameters.jl")
include("parameters_init.jl")
include("parameters_gasc.jl")
include("parameters_actor.jl")

#= ====================== =#
#  PHYSICS and TECHNOLOGY  #
#= ====================== =#
include("physics.jl")
include("technology.jl")

#= ====== =#
#  DDINIT  #
#= ====== =#
include(joinpath("ddinit", "init.jl"))
include(joinpath("ddinit", "init_equilibrium.jl"))
include(joinpath("ddinit", "init_build.jl"))
include(joinpath("ddinit", "init_core_profiles.jl"))
include(joinpath("ddinit", "init_core_sources.jl"))
include(joinpath("ddinit", "init_pf_active.jl"))
include(joinpath("ddinit", "init_others.jl"))

#= ====== =#
#  ACTORS  #
#= ====== =#
abstract type AbstractActor end
function finalize(actor::AbstractActor)
    actor
end
include(joinpath("actors", "equilibrium_actor.jl"))
include(joinpath("actors", "pf_active_actor.jl"))
include(joinpath("actors", "build_actor.jl"))
include(joinpath("actors", "current_actor.jl"))
include(joinpath("actors", "sources_actor.jl"))
include(joinpath("actors", "transport_actor.jl"))

#= ========= =#
#  WORKFLOWS  #
#= ========= =#
include(joinpath("workflows" ,"DB5_validation_workflow.jl"))

#= ====== =#
#= EXPORT =#
#= ====== =#
export step, finalize
export IMAS, evalmissing, @ddtime, constants, ±
export Parameters

end
