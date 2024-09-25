#= ============ =#
#  ActorCurrent  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCurrent{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:SteadyStateCurrent, :QED, :none], "-", "Current actor to run"; default=:SteadyStateCurrent)
    allow_floating_plasma_current::Entry{Bool} = Entry{Bool}("-", "Zero loop voltage if non-inductive fraction exceeds 100% of the target Ip"; default=true)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vloop_from::Switch{Symbol} = switch_get_from(:vloop)
end

mutable struct ActorCurrent{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCurrent{P}
    jt_actor::Union{ActorSteadyStateCurrent{D,P},ActorQED{D,P},ActorNoOperation{D,P}}
end

"""
    ActorCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple ohmic current evolution actors

  - Sets the `j_ohmic`, `j_tor`, `j_total` under `dd.core_profiles.profiles_1d[]`
  - Sets `j_parallel` in `dd.equilibrium.time_slice[].profiles_1d`
  - Updates `bootstrap` and `ohmic` parallel current and heating sources in `dd.core_sources`
"""
function ActorCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCurrent(dd, act.ActorCurrent, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCurrent(dd::IMAS.dd, par::FUSEparameters__ActorCurrent, act::ParametersAllActors; kw...)
    logging_actor_init(ActorCurrent)
    par = par(kw...)
    if par.model == :none
        jt_actor = ActorNoOperation(dd, act.ActorNoOperation)
    elseif par.model == :SteadyStateCurrent
        jt_actor = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent; par.ip_from, par.allow_floating_plasma_current)
    elseif par.model == :QED
        jt_actor = ActorQED(dd, act.ActorQED; par.ip_from, par.vloop_from, par.allow_floating_plasma_current)
    end
    return ActorCurrent(dd, par, jt_actor)
end

"""
    _step(actor::ActorCurrent)

Steps the selected current evolution actor
"""
function _step(actor::ActorCurrent)
    dd = actor.dd

    # freeze jboot and j_non_inductive before updating johmic
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

    # freeze cp1d j_total and j_tor after johmic update
    # important to freeze first j_total and then j_tor
    cp1d = dd.core_profiles.profiles_1d[]
    for field in (:j_total, :j_tor)
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