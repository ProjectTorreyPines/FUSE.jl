#= ============ =#
#  ActorCurrent  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCurrent{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch{Symbol}([:SteadyStateCurrent, :QED], "-", "Current actor to run"; default=:SteadyStateCurrent)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vloop_from::Switch{Symbol} = switch_get_from(:vlopp)
end

mutable struct ActorCurrent{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCurrent{P}
    jt_actor::Union{ActorSteadyStateCurrent{D,P},ActorQED{D,P}}
end

"""
    ActorCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple current evolution actors
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
    if par.model == :SteadyStateCurrent
        jt_actor = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent; par.ip_from)
    elseif par.model == :QED
        jt_actor = ActorQED(dd, act.ActorQED; par.ip_from, par.vloop_from)
    end
    return ActorCurrent(dd, par, jt_actor)
end

"""
    _step(actor::ActorCurrent)

Runs through the selected current evolution actor step
"""
function _step(actor::ActorCurrent)
    step(actor.jt_actor)
    return actor
end

"""
    _finalize(actor::ActorCurrent)

Finalizes the selected current evolution actor finalize
"""
function _finalize(actor::ActorCurrent)
    finalize(actor.jt_actor)
    return actor
end