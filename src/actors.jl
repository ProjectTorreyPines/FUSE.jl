#= =================== =#
#  init IMAS structures #
#= =================== =#
# dummy `init`, `step`, `finalize` simply let IDSs go through them
"""
    dummy `init`
"""
function init(ids::IMAS.IDS)
    return ids
end

"""
    dummy `step`
"""
function step(ids::IMAS.IDS)
    return ids
end

"""
    dummy `finalize`
"""
function finalize(ids::IMAS.IDS)
    return ids
end

#= ============ =#
#  AbstractActor #
#= ============ =#
abstract type AbstractActor end

"""
    Take a step with a given actor
"""
function Base.step(actor::AbstractActor)
    error("Function step() not defined for actor of type $(typeof(actor))")
end

"""
    store output data in IDS
    NOTE: actors should take dd and output dd, operating in place
"""
function finalize(actor::AbstractActor)
    error("Function finalize() not defined for actor of type $(typeof(actor))")
end

#= =========== =#
#  Equilibrium  #
#= =========== =#

abstract type EquilibriumActor <: AbstractActor end

include("actors/equilibrium_actor.jl")

export SolovevEquilibriumActor

#= ===== =#
#  Coils  #
#= ===== =#

abstract type CoilsActor <: AbstractActor end

include("actors/coils_actor.jl")

export PFcoilsOptActor

#= ============ =#
#  Radial Build  #
#= ============ =#

abstract type RadialBuildActor <: AbstractActor end

include("actors/radial_build_actor.jl")

