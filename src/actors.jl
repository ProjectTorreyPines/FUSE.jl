#= =================== =#
#  init IMAS structures #
#= =================== =#

function init(ids::IMAS.IDS, time::Real)
    error("Function init() not defined for ids of type $(typeof(ids))")
end

function init(equilibrium::IMAS.equilibrium, time::Real=0.0;
        B0, R0, ϵ, δ, κ, beta_n, ip,
        diverted::Bool=false,
        xpoint::Union{NTuple{2},Nothing}=(diverted ? (R0 * (1 - 1.1 * δ * ϵ), -R0 * 1.1 * κ * ϵ) : nothing),
        symmetric::Bool=(xpoint === nothing))
    # all this should be switched to use 'boundary'
    time_index = get_time_index(equilibrium.time_slice, time)
    eqt = equilibrium.time_slice[time_index]
    eqt.boundary.minor_radius = ϵ * R0
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.elongation = κ
    eqt.boundary.triangularity = δ
    set_field_time_array(equilibrium.vacuum_toroidal_field, :b0, time_index, B0)
    equilibrium.vacuum_toroidal_field.r0 = R0
    eqt.global_quantities.ip = ip
    eqt.global_quantities.beta_normal = beta_n
    return equilibrium
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

include("actors/solovev_equilibrium_actor.jl")

