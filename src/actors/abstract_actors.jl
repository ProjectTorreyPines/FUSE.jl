abstract type AbstractActor end
abstract type FacilityAbstractActor <: AbstractActor end
abstract type ReactorAbstractActor <: AbstractActor end
abstract type HCDAbstractActor <: AbstractActor end
abstract type PlasmaAbstractActor <: AbstractActor end

function _finalize(actor::AbstractActor)
    actor
end

function (actor::AbstractActor)(;kwargs...)
    step(actor)
    finalize(actor)
end