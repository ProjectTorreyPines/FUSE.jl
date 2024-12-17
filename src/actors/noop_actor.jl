#= ================ =#
#  ActorNoOperation  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorNoOperation{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
end

mutable struct ActorNoOperation{D,P} <: AbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorNoOperation{P}
    function ActorNoOperation(dd::IMAS.dd{D}, par::FUSEparameters__ActorNoOperation{P}; kw...) where {D<:Real,P<:Real}
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorNoOperation(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor doesn't do anything. It can be useful to turn off models, for example.
"""
function ActorNoOperation(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorNoOperation(dd, act.ActorNoOperation; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorNoOperation)
    return actor
end

function _finalize(actor::ActorNoOperation)
    return actor
end
