#= =================== =#
#  ActorCoreTransport  #
#= =================== =#
mutable struct ActorCoreTransport <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    turb_actor::PlasmaAbstractActor
    neoclassical_actor::PlasmaAbstractActor
end

Base.@kwdef mutable struct FUSEparameters__ActorCoreTransport{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    rho_transport = Entry(AbstractVector, "", "rho core transport grid"; default=0.2:0.1:0.8)
    turbulence_actor = Switch(Symbol, [:TGLF, :None], "", "Turbulence Actor to run"; default=:TGLF)
    neoclassical_actor = Switch(Symbol, [:Neoclassical, :None], "", "Neocalssical actor to run"; default=:Neoclassical)
end

"""
    ActorCoreTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple equilibrium actors
"""
function ActorCoreTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCoreTransport(kw...)
    actor = ActorCoreTransport(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCoreTransport(dd::IMAS.dd, par::ParametersActor, act::ParametersAllActors; kw...)
    par = par(kw...)

    if par.turbulence_actor == :TGLF
        act.ActorTGLF.rho_transport = par.rho_transport
        turb_actor = ActorTGLF(dd, act.ActorTGLF)
    else
        error("turbulence_actor $(par.turbulence_actor) is not supported yet")
    end

    if par.neoclassical_actor == :Neoclassical
        act.ActorNeoclassical.rho_transport = par.rho_transport
        neoclassical_actor = ActorNeoclassical(dd, act.ActorNeoclassical)
    else
        error("neoclassical_actor $(par.neoclassical_actor) is not supported yet")
    end

    return ActorCoreTransport(dd, par, turb_actor, neoclassical_actor)
end

"""
    step(actor::ActorCoreTransport)

Runs through the selected equilibrium actor's step
"""
function _step(actor::ActorCoreTransport)
    step(actor.turb_actor)
    step(actor.neoclassical_actor)
    return actor
end

"""
    finalize(actor::ActorCoreTransport)

Finalizes the selected equilibrium actor
"""
function finalize(actor::ActorCoreTransport)
    finalize(actor.turb_actor)
    finalize(actor.neoclassical_actor)
    return actor
end