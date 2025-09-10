#= ======== =#
#  ActorSOL  #
#= ======== =#
@actor_parameters_struct ActorSOL{T} begin
    model::Switch{Symbol} = Switch{Symbol}([:box, :replay, :none], "-", "SOL actor to run"; default=:box)
end

mutable struct ActorSOL{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSOL{P}}
    act::ParametersAllActors{P}
    SOL_actor::Union{Nothing,ActorSOLBox{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
end

"""
    ActorSOL(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run different SOL actors
"""
function ActorSOL(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSOL(dd, act.ActorSOL, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorSOL(dd::IMAS.dd, par::FUSEparameters__ActorSOL, act::ParametersAllActors; kw...)
    logging_actor_init(ActorSOL)
    par = OverrideParameters(par; kw...)

    if par.model == :box
        SOL_actor = ActorSOLBox(dd, act.ActorSOLBox)
    elseif par.model == :replay
        SOL_actor = ActorReplay(dd, act.ActorReplay, actor)
    elseif par.model == :none
        SOL_actor = ActorNoOperation(dd, act.ActorNoOperation)
    end

    return ActorSOL(dd, par, act, SOL_actor)
end

function _step(actor::ActorSOL)
    step(actor.SOL_actor)
    return actor
end

function _finalize(actor::ActorSOL)
    finalize(actor.SOL_actor)
    return actor
end
