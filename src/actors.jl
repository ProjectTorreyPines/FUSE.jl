#= =================== =#
#  init IMAS structures #
#= =================== =#
"""
    dummy `init`
"""
function init(ids::IMAS.IDS)
    error("Function init() not defined for IDS of type $(typeof(ids))")
end


#= ============ =#
#  AbstractActor #
#= ============ =#
abstract type AbstractActor end

"""
    Base.step(actor::AbstractActor)

Placeholder function to take a step with a given actor
"""
function Base.step(actor::AbstractActor)
    error("Function step() not defined for actor of type $(typeof(actor))")
end

"""
    finalize(actor::AbstractActor)

Dummy function to finalize actor and store its in a IDS
"""
function finalize(actor::AbstractActor)
    actor
end

#= =========== =#
#  Equilibrium  #
#= =========== =#
include("actors/equilibrium_actor.jl")

export SolovevEquilibriumActor

#= ===== =#
#  Coils  #
#= ===== =#
include("actors/coils_actor.jl")

export PFcoilsOptActor

#= ============ =#
#  Radial Build  #
#= ============ =#
include("actors/radial_build_actor.jl")

#= ========= =#
#  Transport  #
#= ========= =#

include("actors/transport_actor.jl")
