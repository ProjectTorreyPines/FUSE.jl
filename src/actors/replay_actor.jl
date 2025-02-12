#= =========== =#
#  ActorReplay  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorReplay{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    replay_dd::Entry{IMAS.dd{T}} = Entry{IMAS.dd{T}}("-", "`dd` to replay data from")
end

mutable struct ActorReplay{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorReplay{P}
    replay_dd::IMAS.dd{D}
    base_actor::AbstractActor{D,P}
end

function ActorReplay(dd::IMAS.dd{D}, par::FUSEparameters__ActorReplay{P}, base_actor::AbstractActor{D,P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorReplay)
    # NOTE: we don't do `par = par(kw...)` here to avoid the deepcopy of an entire `dd`
    if :replay_dd in keys(kw)
        replay_dd = kw.replay_dd
    elseif !ismissing(par, :replay_dd)
        replay_dd = par.replay_dd
    else
        replay_dd = IMAS.dd{D}()
    end
    return ActorReplay{D,P}(dd, par, replay_dd, base_actor)
end

"""
    ActorReplay(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor is meant to be used as a specific model in generic actors
for replaying past behaviour from data specified in `act.ActorReplay.replay_dd`.

The behaviour of this actor (what it replays) is defined by the dispatch on these functions:

  - `_step(::ActorReplay, actor::Actor???, replay_dd::IMAS.dd)`
  - `_finalize(::ActorReplay, actor::Actor???, replay_dd::IMAS.dd)`
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
