#= =================== =#
#  ActorFluxCalculator  #
#= =================== =#
@actor_parameters_struct ActorFluxCalculator{T} begin
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho core transport grid"; default=0.25:0.1:0.85)
    turbulence_model::Switch{Symbol} = Switch{Symbol}([:TGLF, :QLGYRO, :analytic, :none], "-", "Turbulence model to use"; default=:TGLF)
    neoclassical_model::Switch{Symbol} = Switch{Symbol}([:neoclassical, :none], "-", "Neocalssical model to use"; default=:neoclassical)
end

mutable struct ActorFluxCalculator{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorFluxCalculator{P}}
    act::ParametersAllActors{P}
    actor_turb::Union{ActorTGLF{D,P},ActorQLGYRO{D,P},ActorAnalyticTurbulence{D,P},ActorNoOperation{D,P}}
    actor_neoc::Union{ActorNeoclassical{D,P},ActorNoOperation{D,P}}
end

"""
    ActorFluxCalculator(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a unified interface to run turbulent and neoclassical transport model actors.

This compound actor manages the execution of both turbulent transport models 
(TGLF, QLGYRO, analytic, or none) and neoclassical transport models (neoclassical or none).
The actor coordinates the calculation of transport fluxes from both physics mechanisms
and ensures they are properly stored in `dd.core_transport`.

Turbulence model options:
- `:TGLF`: TGLF-based models (TGLF, TGLFNN, GKNN, TJLF)
- `:QLGYRO`: Quasi-linear gyrokinetic transport via CGYRO
- `:analytic`: Simple analytic transport models
- `:none`: No turbulent transport

Neoclassical model options:
- `:neoclassical`: Collisional transport (Chang-Hinton, NEO, Hirshman-Sigmar)
- `:none`: No neoclassical transport
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
    elseif par.turbulence_model == :analytic
        actor_turb = ActorAnalyticTurbulence(dd, act.ActorAnalyticTurbulence; par.rho_transport)
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

Runs the selected turbulent and neoclassical transport actors in sequence.

First executes the turbulent transport actor, then the neoclassical transport actor.
Each actor calculates its respective transport fluxes and stores them in the
appropriate sections of `dd.core_transport`.
"""
function _step(actor::ActorFluxCalculator)
    step(actor.actor_turb)
    step(actor.actor_neoc)
    return actor
end

"""
    _finalize(actor::ActorFluxCalculator)

Finalizes both the turbulent and neoclassical transport actors.

Calls the finalize methods of both sub-actors to ensure their results are properly
written to the `dd.core_transport` data structure.
"""
function _finalize(actor::ActorFluxCalculator)
    finalize(actor.actor_turb)
    finalize(actor.actor_neoc)
    return actor
end