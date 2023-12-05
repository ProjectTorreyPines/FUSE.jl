#= =================== =#
#  ActorFluxCalculator  #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxCalculator{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho core transport grid"; default=0.25:0.1:0.85)
    turbulence_model::Switch{Symbol} = Switch{Symbol}([:TGLF, :none], "-", "Turbulence model to use"; default=:TGLF)
    neoclassical_model::Switch{Symbol} = Switch{Symbol}([:neoclassical, :none], "-", "Neocalssical model to use"; default=:neoclassical)
end

mutable struct ActorFluxCalculator{D,P} <: PlasmaAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFluxCalculator{P}
    actor_turb::Union{ActorTGLF{D,P},ActorNoOperation{D,P}}
    actor_neoc::Union{ActorNeoclassical{D,P},ActorNoOperation{D,P}}
end

"""
    ActorFluxCalculator(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple transport model actors
"""
function ActorFluxCalculator(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFluxCalculator(dd, act.ActorFluxCalculator, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFluxCalculator(dd::IMAS.dd, par::FUSEparameters__ActorFluxCalculator, act::ParametersAllActors; kw...)
    par = par(kw...)

    if par.turbulence_model == :none
        logging(Logging.Debug, :actors, "ActorFluxCalculator: turbulent transport disabled")
        actor_turb = ActorNoOperation(dd, act.ActorNoOperation)
    elseif par.turbulence_model == :TGLF
        act.ActorTGLF.rho_transport = par.rho_transport
        actor_turb = ActorTGLF(dd, act.ActorTGLF)
    end

    if par.neoclassical_model == :none
        logging(Logging.Debug, :actors, "ActorFluxCalculator: neoclassical transport disabled")
        actor_neoc = ActorNoOperation(dd, act.ActorNoOperation)
    elseif par.neoclassical_model == :neoclassical
        act.ActorNeoclassical.rho_transport = par.rho_transport
        actor_neoc = ActorNeoclassical(dd, act.ActorNeoclassical)
    end

    return ActorFluxCalculator(dd, par, actor_turb, actor_neoc)
end

"""
    _step(actor::ActorFluxCalculator)

Runs through the selected equilibrium actor's step
"""
function _step(actor::ActorFluxCalculator)
    step(actor.actor_turb)
    step(actor.actor_neoc)
    return actor
end

"""
    _finalize(actor::ActorFluxCalculator)

Finalizes the selected equilibrium actor
"""
function _finalize(actor::ActorFluxCalculator)
    finalize(actor.actor_turb)
    finalize(actor.actor_neoc)
    return actor
end