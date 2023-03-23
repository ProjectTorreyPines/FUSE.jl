abstract type AbstractActor end
abstract type FacilityAbstractActor <: AbstractActor end
abstract type ReactorAbstractActor <: AbstractActor end
abstract type HCDAbstractActor <: AbstractActor end
abstract type PlasmaAbstractActor <: AbstractActor end

function logging_actor_init(typeof_actor::DataType, args...; kw...)
    logging(Logging.Debug, :actors, "$typeof_actor @ init")
end

#= ==== =#
#  step  #
#= ==== =#
"""
    step(actor::T, args...; kw...) where {T<:AbstractActor}

Run a actor
"""
function step(actor::T, args...; kw...) where {T<:AbstractActor}
    logging(Logging.Info, :actors, "$(typeof(actor)) @ step")
    timer_name = replace(string(typeof(actor).name.name),r"^Actor" => "")
    TimerOutputs.reset_timer!(timer_name)
    TimerOutputs.@timeit timer timer_name begin
        _step(actor, args...; kw...)::T
    end
    return actor
end

#= ======== =#
#  finalize  #
#= ======== =#
function _finalize(actor::AbstractActor)
    return actor
end

"""
    finalize(actor::T, args...; kw...) where {T<:AbstractActor}

Finalize the actor run. This is typically used to update dd.
"""
function finalize(actor::T, args...; kw...) where {T<:AbstractActor}
    logging(Logging.Debug, :actors, "$(typeof(actor)) @finalize")
    _finalize(actor, args...; kw...)::T
    return actor
end

#= ======= =#
#  prepare  #
#= ======= =#
"""
    prepare(actor_type::DataType, dd::IMAS.dd, act::ParametersAllActors; kw...)

Dispatch `prepare` function for different actors based on actor_type that is passed
"""
function prepare(dd::IMAS.dd, actor_name::Symbol, act::ParametersAllActors; kw...)
    prepare(dd, Val{actor_name}, act; kw...)
    return dd
end