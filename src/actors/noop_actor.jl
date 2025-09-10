#= ================ =#
#  ActorNoOperation  #
#= ================ =#
@actor_parameters_struct ActorNoOperation{T} begin
end

mutable struct ActorNoOperation{D,P} <: AbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorNoOperation{P}}
    function ActorNoOperation(dd::IMAS.dd{D}, par::FUSEparameters__ActorNoOperation{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorNoOperation)
        par = OverrideParameters(par; kw...)
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
