#= =========== =#
#  ActorReplay  #
#= =========== =#
@actor_parameters_struct ActorReplay{T} begin
    replay_dd::Entry{IMAS.DD} = Entry{IMAS.DD}("-", "`dd` to replay data from")
end

mutable struct ActorReplay{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorReplay{P}}
    replay_dd::IMAS.dd{D}
    base_actor::AbstractActor{D,P}
end

function ActorReplay(dd::IMAS.dd{D}, par::FUSEparameters__ActorReplay{P}, base_actor::AbstractActor{D,P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorReplay)
    par = OverrideParameters(par; kw...)
    if ismissing(par, :replay_dd)
        replay_dd = IMAS.dd{D}()
    elseif typeof(par.replay_dd) == typeof(dd)
        replay_dd = par.replay_dd
    else
        replay_dd = IMAS.dd{D}()
        fill!(replay_dd, par.replay_dd)
    end
    return ActorReplay{D,P}(dd, par, replay_dd, base_actor)
end

function ActorReplay(dd::IMAS.dd{D}, par::FUSEparameters__ActorReplay{P}; kw...) where {D<:Real,P<:Real}
    # This constructor is only used to conform with the standard SingleAbstractActor call signature
    # which is necessary to not break the generation of the FUSE documentation.
    ActorReplay(dd, par, ActorNoOperation(dd, ParametersAllActors{D,P}()))
end

"""
    ActorReplay(dd::IMAS.dd, act::ParametersAllActors; kw...)

Replays plasma profiles and behavior from experimental data or previous simulation results.

This actor serves as a transport model that imposes predetermined plasma profiles
instead of calculating them from physics models. It reads data from `act.ActorReplay.replay_dd`
and applies it to the current simulation state.

The specific replay behavior depends on the context where this actor is used, defined by
specialized dispatch methods:
- `_step(::ActorReplay, actor::Actor???, replay_dd::IMAS.dd)` - defines what data to replay
- `_finalize(::ActorReplay, actor::Actor???, replay_dd::IMAS.dd)` - handles post-processing

Common use cases:
- Replaying experimental temperature/density profiles in flux matching
- Using previous simulation results as boundary conditions  
- Imposing prescribed time evolution for sensitivity studies
- Validation against experimental data

The replay data source is specified via the `replay_dd` parameter, which can be
an experimental shot, previous simulation, or synthetic data.
"""
function ActorReplay(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorReplay(dd, act.ActorReplay; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorReplay)
    return _step(actor, actor.base_actor, actor.replay_dd)
end

function _finalize(actor::ActorReplay)
    return _finalize(actor, actor.base_actor, actor.replay_dd)
end

function _step(replay_actor::ActorReplay, actor::AbstractActor, replay_dd::IMAS.dd)
    return replay_actor
end

function _finalize(replay_actor::ActorReplay, actor::AbstractActor, replay_dd::IMAS.dd)
    return replay_actor
end
