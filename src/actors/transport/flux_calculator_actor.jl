#= =================== =#
#  ActorFluxCalculator  #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxCalculator{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho core transport grid"; default=0.25:0.1:0.85)
    turbulence_model::Switch{Symbol} = Switch{Symbol}([:TGLF, :QLGYRO, :none], "-", "Turbulence model to use"; default=:TGLF)
    neoclassical_model::Switch{Symbol} = Switch{Symbol}([:neoclassical, :none], "-", "Neocalssical model to use"; default=:neoclassical)
end

mutable struct ActorFluxCalculator{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorFluxCalculator{P}}
    act::ParametersAllActors{P}
    actor_turb::Union{ActorTGLF{D,P},ActorQLGYRO{D,P},ActorAnalyticalTurbulence{D,P},ActorNoOperation{D,P}}
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
    logging_actor_init(ActorFluxCalculator)
    par = OverrideParameters(par; kw...)

    if par.turbulence_model == :none
        actor_turb = ActorNoOperation(dd, act.ActorNoOperation)
    elseif par.turbulence_model == :TGLF
        actor_turb = ActorTGLF(dd, act.ActorTGLF; par.rho_transport)
    elseif par.turbulence_model == :QLGYRO
        actor_turb = ActorQLGYRO(dd, act.ActorQLGYRO; par.rho_transport)
    elseif par.turbulence_model == :Analytical
        actor_turb = ActorAnalyticalTurbulence(dd, act.ActorAnalyticalTurbulence; par.rho_transport)
    end

    if par.neoclassical_model == :none
        actor_neoc = ActorNoOperation(dd, act.ActorNoOperation)
    elseif par.neoclassical_model == :neoclassical
        act.ActorNeoclassical.rho_transport = par.rho_transport
        actor_neoc = ActorNeoclassical(dd, act.ActorNeoclassical)
    end

    return ActorFluxCalculator(dd, par, act, actor_turb, actor_neoc)
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