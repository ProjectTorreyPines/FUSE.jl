#= ============ =#
#  ActorCurrent  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCurrent{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:SteadyStateCurrent, :QED, :replay, :none], "-", "Current actor to run"; default=:QED)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vloop_from::Switch{Symbol} = switch_get_from(:vloop)
end

mutable struct ActorCurrent{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorCurrent{P}}
    act::ParametersAllActors{P}
    jt_actor::Union{ActorSteadyStateCurrent{D,P},ActorQED{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
end

"""
    ActorCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple ohmic current evolution actors

  - The fundamental quantitiy being solved is `j_total` in `dd.core_profiles.profiles_1d[]`
  - Also updates:

      + `j_ohmic`, `j_tor` in `dd.core_profiles.profiles_1d[]`
      + `j_parallel` in `dd.equilibrium.time_slice[].profiles_1d`
      + `bootstrap` and `ohmic` parallel current and heating sources in `dd.core_sources`
"""
function ActorCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCurrent(dd, act.ActorCurrent, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCurrent(dd::IMAS.dd, par::FUSEparameters__ActorCurrent, act::ParametersAllActors; kw...)
    logging_actor_init(ActorCurrent)
    par = OverrideParameters(par; kw...)

    noop = ActorNoOperation(dd, act.ActorNoOperation)
    actor = ActorCurrent(dd, par, act, noop)

    if par.model == :SteadyStateCurrent
        actor.jt_actor = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent; par.ip_from)
    elseif par.model == :QED
        actor.jt_actor = ActorQED(dd, act.ActorQED, act; par.ip_from, par.vloop_from)
    elseif par.model == :replay
        actor.jt_actor = ActorReplay(dd, act.ActorReplay, actor)
    end

    return actor
end

"""
    _step(actor::ActorCurrent)

Steps the selected current evolution actor
"""
function _step(actor::ActorCurrent)
    dd = actor.dd

    # freeze jboot and j_non_inductive before updating j_ohmic
    cp1d = dd.core_profiles.profiles_1d[]
    for field in [:j_bootstrap, :j_non_inductive]
        IMAS.refreeze!(cp1d, field)
    end

    step(actor.jt_actor)

    return actor
end

"""
    _finalize(actor::ActorCurrent)

Finalizes the selected current evolution actor
"""
function _finalize(actor::ActorCurrent)
    dd = actor.dd

    finalize(actor.jt_actor)

    # freeze cp1d j_ohmic and j_tor after j_total update
    cp1d = dd.core_profiles.profiles_1d[]
    for field in (:j_ohmic, :j_tor)
        IMAS.refreeze!(cp1d, field)
    end

    # similarly, freeze parallel plasma current
    eqt = dd.equilibrium.time_slice[]
    IMAS.refreeze!(eqt.profiles_1d, :j_parallel)

    # update core_sources related to current
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)

    return actor
end

function _step(replay_actor::ActorReplay, actor::ActorCurrent, replay_dd::IMAS.dd)
    dd = actor.dd

    time0 = dd.global_time
    cp1d = dd.core_profiles.profiles_1d[time0]
    replay_cp1d = replay_dd.core_profiles.profiles_1d[time0]

    cp1d.j_total = replay_cp1d.j_total

    return replay_actor
end
