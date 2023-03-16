abstract type AbstractActor end
abstract type FacilityAbstractActor <: AbstractActor end
abstract type ReactorAbstractActor <: AbstractActor end
abstract type HCDAbstractActor <: AbstractActor end
abstract type PlasmaAbstractActor <: AbstractActor end

function logging_actor_init(typeof_actor::DataType, args...; kw...)
    logging(Logging.Debug, :actors, "$typeof_actor @ init")
end

function step(actor::T, args...; kw...) where {T<:AbstractActor}
    logging(Logging.Info, :actors, "$(typeof(actor)) @ step")
    timer_name = replace(string(typeof(actor).name.name),r"^Actor" => "")
    TimerOutputs.reset_timer!(to, timer_name)
    TimerOutputs.@timeit to timer_name begin
        _step(actor, args...; kw...)::T
    end
    return actor
end

function _finalize(actor::AbstractActor)
    return actor
end

function finalize(actor::T, args...; kw...) where {T<:AbstractActor}
    logging(Logging.Debug, :actors, "$(typeof(actor)) @finalize")
    _finalize(actor, args...; kw...)::T
    return actor
end