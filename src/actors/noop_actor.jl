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

A no-operation actor that performs no calculations or modifications.

This actor serves as a placeholder when specific physics models need to be disabled
in compound actors or actor workflows. It implements the standard actor interface
(step and finalize methods) but performs no operations, making it useful for:

- Disabling transport models in ActorFluxCalculator
- Turning off specific physics in complex actor chains  
- Serving as a default/null option in actor selection logic
- Testing and debugging actor workflows

The actor simply returns itself from both `_step()` and `_finalize()` methods.
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
