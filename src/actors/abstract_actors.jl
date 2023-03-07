abstract type AbstractActor end
abstract type FacilityAbstractActor <: AbstractActor end
abstract type ReactorAbstractActor <: AbstractActor end
abstract type HCDAbstractActor <: AbstractActor end
abstract type PlasmaAbstractActor <: AbstractActor end

function logging_actor_init(typeof_actor::DataType, args...; kw...)
    logging(Logging.Debug, :actors, "$typeof_actor @ init")
end

function step(actor::AbstractActor, args...; kw...)
    logging(Logging.Info, :actors, "$(typeof(actor)) @ step")
    return _step(actor, args...; kw...)
end

function _finalize(actor::AbstractActor)
    actor
end

function finalize(actor::AbstractActor, args...; kw...)
    logging(Logging.Debug, :actors, "$(typeof(actor)) @finalize")
    return _finalize(actor, args...; kw...)
end