#= =================== =#
#  ActorFluxCalculator  #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxCalculator{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho core transport grid"; default=0.2:0.1:0.8)
    turbulence_model::Switch{Symbol} = Switch{Symbol}([:TGLF, :none], "-", "Turbulence model to use"; default=:TGLF)
    neoclassical_model::Switch{Symbol} = Switch{Symbol}([:neoclassical, :none], "-", "Neocalssical model to use"; default=:neoclassical)
end

mutable struct ActorFluxCalculator{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFluxCalculator{P}
    turb_actor::ActorTGLF{D,P}
    neoc_actor::ActorNeoclassical{D,P}
end

"""
    ActorFluxCalculator(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple equilibrium actors
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
    elseif par.turbulence_model == :TGLF
        act.ActorTGLF.rho_transport = par.rho_transport
        turb_actor = ActorTGLF(dd, act.ActorTGLF)
    end

    if par.neoclassical_model == :none
        logging(Logging.Debug, :actors, "ActorFluxCalculator: neoclassical transport disabled")
    elseif par.neoclassical_model == :neoclassical
        act.ActorNeoclassical.rho_transport = par.rho_transport
        neoc_actor = ActorNeoclassical(dd, act.ActorNeoclassical)
    end

    return ActorFluxCalculator(dd, par, turb_actor, neoc_actor)
end

"""
    _step(actor::ActorFluxCalculator)

Runs through the selected equilibrium actor's step
"""
function _step(actor::ActorFluxCalculator)
    step(actor.turb_actor)
    step(actor.neoc_actor)
    return actor
end

"""
    _finalize(actor::ActorFluxCalculator)

Finalizes the selected equilibrium actor
"""
function _finalize(actor::ActorFluxCalculator)
    finalize(actor.turb_actor)
    finalize(actor.neoc_actor)
    return actor
end