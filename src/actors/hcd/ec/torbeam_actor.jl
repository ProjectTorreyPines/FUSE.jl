import TORBEAM
using Plots

#= =========== =#
#  ActorTORBEAM  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorTORBEAM{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
end

mutable struct ActorTORBEAM{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTORBEAM{P}
    torbeam_params::TORBEAM.TorbeamParams
end

function ActorTORBEAM(dd::IMAS.dd, par::FUSEparameters__ActorTORBEAM; kw...)
    logging_actor_init(ActorTORBEAM)
    par = par(kw...)
    return ActorTORBEAM(dd, par, TORBEAM.TorbeamParams())
end

"""
    ActorTORBEAM(dd::IMAS.dd, act::ParametersAllActors; kw...)
"""
function ActorTORBEAM(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTORBEAM(dd, act.ActorTORBEAM; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorTORBEAM)
    dd = actor.dd
    TORBEAM.torbeam!(dd, actor.torbeam_params)
    return actor
end
