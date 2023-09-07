#= ========= =#
#  ActorNoOP  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorNoOP{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorNoOP{D,P} <: AbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorNoOP{P}
    function ActorNoOP(dd::IMAS.dd{D}, par::FUSEparameters__ActorNoOP{P}; kw...) where {D<:Real,P<:Real}
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorNoOP(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor doesn't do anything. It can be useful to turn off models, for example.
"""
function ActorNoOP(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorNoOP(dd, act.ActorNoOP; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorNoOP)
    return actor
end

function _finalize(actor::ActorNoOP)
    return actor
end

function step(actor::ActorNoOP)
    return actor
end

function finalize(actor::ActorNoOP)
    return actor
end