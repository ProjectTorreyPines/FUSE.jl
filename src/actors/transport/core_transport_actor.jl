#= =================== =#
#  ActorCoreTransport  #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorCoreTransport{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    rho_transport::Entry{AbstractVector{<:T}} = Entry(AbstractVector{<:T}, "-", "rho core transport grid"; default=0.2:0.1:0.8)
    turbulence_model::Switch{Symbol} = Switch(Symbol, [:TGLF, :none], "-", "Turbulence model to use"; default=:TGLF)
    neoclassical_model::Switch{Symbol} = Switch(Symbol, [:neoclassical, :none], "-", "Neocalssical model to use"; default=:neoclassical)
end

mutable struct ActorCoreTransport <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorCoreTransport
    turb_actor::PlasmaAbstractActor
    neoc_actor::PlasmaAbstractActor
end

"""
    ActorCoreTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple equilibrium actors
"""
function ActorCoreTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCoreTransport
    actor = ActorCoreTransport(dd, par, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCoreTransport(dd::IMAS.dd, par::FUSEparameters__ActorCoreTransport, act::ParametersAllActors; kw...)
    par = par(kw...)

    if par.turbulence_model == :none
        logging(Logging.Debug, :actors, "ActorCoreTransport: turbulent transport disabled")
    elseif par.turbulence_model == :TGLF
        act.ActorTGLF.rho_transport = par.rho_transport
        turb_actor = ActorTGLF(dd, act.ActorTGLF)
    else
        error("turbulence_model `$(par.turbulence_model)` is not supported yet")
    end

    if par.neoclassical_model == :none
        logging(Logging.Debug, :actors, "ActorCoreTransport: neoclassical transport disabled")
    elseif par.neoclassical_model == :neoclassical
        act.ActorNeoclassical.rho_transport = par.rho_transport
        neoc_actor = ActorNeoclassical(dd, act.ActorNeoclassical)
    else
        error("neoclassical_model `$(par.neoclassical_model)` is not supported yet")
    end

    return ActorCoreTransport(dd, par, turb_actor, neoc_actor)
end

"""
    step(actor::ActorCoreTransport)

Runs through the selected equilibrium actor's step
"""
function _step(actor::ActorCoreTransport)
    step(actor.turb_actor)
    step(actor.neoc_actor)
    return actor
end

"""
    finalize(actor::ActorCoreTransport)

Finalizes the selected equilibrium actor
"""
function _finalize(actor::ActorCoreTransport)
    finalize(actor.turb_actor)
    finalize(actor.neoc_actor)
    return actor
end