#= === =#
#  NBI  #
#= === =#
Base.@kwdef mutable struct FUSEparameters__ActorNBI{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:simple, :RABBIT], "-", "NBI model"; default=:simple)
end

mutable struct ActorNBI{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorNBI{P}
    act::ParametersAllActors
    nbi_actor::Union{ActorSimpleNB{D,P},ActorRABBIT{D,P}}
end

"""
    ActorNBI(dd::IMAS.dd, act::ParametersAllActors; kw...)

"""
function ActorNBI(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorNBI(dd, act.ActorNBI, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorNBI(dd::IMAS.dd, par::FUSEparameters__ActorNBI, act::ParametersAllActors; kw...)
    logging_actor_init(ActorNBI)

    par = par(kw...)

    if par.model == :simple 
        nbi_actor = ActorSimpleNB(dd, act.ActorSimpleNB)
    elseif par.model == :RABBIT
        nbi_actor = ActorRABBIT(dd, act.ActorRABBIT)
    end

    return ActorNBI(dd, par, act, nbi_actor)
end

function _step(actor::ActorNBI)
    step(actor.nbi_actor)
    return actor
end

function _finalize(actor::ActorNBI)
    finalize(actor.nbi_actor)
    return actor 
end