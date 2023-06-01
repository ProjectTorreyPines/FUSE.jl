#= ======== =#
#  ActorHCD  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorHCD{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    ec_model::Switch{Symbol} = Switch{Symbol}([:ECsimple], "-", "EC source actor to run"; default=:ECsimple)
    ic_model::Switch{Symbol} = Switch{Symbol}([:ICsimple], "-", "IC source actor to run"; default=:ICsimple)
    lh_model::Switch{Symbol} = Switch{Symbol}([:LHsimple], "-", "LH source actor to run"; default=:LHsimple)
    nb_model::Switch{Symbol} = Switch{Symbol}([:NBsimple], "-", "NB source actor to run"; default=:NBsimple)
end

mutable struct ActorHCD{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorHCD{P}
    ec_actor::ActorECsimple{D,P}
    ic_actor::ActorICsimple{D,P}
    lh_actor::ActorLHsimple{D,P}
    nb_actor::ActorNBsimple{D,P}
end

"""
    ActorHCD(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run HCD actors
"""
function ActorHCD(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorHCD(dd, act.ActorHCD, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorHCD(dd::IMAS.dd, par::FUSEparameters__ActorHCD, act::ParametersAllActors; kw...)
    logging_actor_init(ActorHCD)
    par = par(kw...)
    if par.ec_model == :ECsimple
        ec_actor = ActorECsimple(dd, act.ActorECsimple)
    end
    if par.ic_model == :ICsimple
        ic_actor = ActorICsimple(dd, act.ActorICsimple)
    end
    if par.lh_model == :LHsimple
        lh_actor = ActorLHsimple(dd, act.ActorLHsimple)
    end
    if par.nb_model == :NBsimple
        nb_actor = ActorNBsimple(dd, act.ActorNBsimple)
    end
    return ActorHCD(dd, par, ec_actor, ic_actor, lh_actor, nb_actor)
end

"""
    _step(actor::ActorHCD)

Runs through the selected HCD actor's step
"""
function _step(actor::ActorHCD)
    step(actor.ec_actor)
    step(actor.ic_actor)
    step(actor.lh_actor)
    step(actor.nb_actor)
    return actor
end

"""
    finalize(actor::ActorHCD)

Finalizes the selected CHD actor's finalize
"""
function _finalize(actor::ActorHCD)
    finalize(actor.ec_actor)
    finalize(actor.ic_actor)
    finalize(actor.lh_actor)
    finalize(actor.nb_actor)
    return actor
end