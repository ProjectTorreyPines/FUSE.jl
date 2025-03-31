import GACODE

#= ======================= =#
#  ActorAnalyticTurbulence  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorAnalyticTurbulence{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:GyroBohm], "-", "Analytic transport model"; default=:GyroBohm)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute fluxes on"; default=0.25:0.1:0.85)
end

mutable struct ActorAnalyticTurbulence{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorAnalyticTurbulence{P}}
    flux_solutions::Vector{<:GACODE.FluxSolution}
end

"""
    ActorAnalyticTurbulence(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates analytic turbulence models
"""
function ActorAnalyticTurbulence(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorAnalyticTurbulence(dd, act.ActorAnalyticTurbulence; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorAnalyticTurbulence(dd::IMAS.dd, par::FUSEparameters__ActorAnalyticTurbulence; kw...)
    logging_actor_init(ActorAnalyticTurbulence)
    par = OverrideParameters(par; kw...)
    return ActorAnalyticTurbulence(dd, par, GACODE.FluxSolution[])
end

"""
    _step(actor::ActorAnalyticTurbulence)

Runs analytic turbulent transport model on a vector of gridpoints
"""
function _step(actor::ActorAnalyticTurbulence)
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]

    flux_solution = GACODE.FluxSolution(1.0, 1.0, 1.0, [1.0 for ion in cp1d.ion], 1.0)
    actor.flux_solutions = [flux_solution for irho in par.rho_transport]

    return actor
end

"""
    _finalize(actor::ActorAnalyticTurbulence)

Writes results to dd.core_transport
"""
function _finalize(actor::ActorAnalyticTurbulence)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = string(par.model)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux,  :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end
