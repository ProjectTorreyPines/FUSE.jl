import GACODE

#= ========================= =#
#  ActorAnalyticalTurbulence  #
#= ========================= =#
Base.@kwdef mutable struct FUSEparameters__ActorAnalyticalTurbulence{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:GyroBohm], "-", "Analytical transport model"; default=:GyroBohm)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute fluxes on"; default=0.25:0.1:0.85)
end

mutable struct ActorAnalyticalTurbulence{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorAnalyticalTurbulence{P}}
    flux_solutions::Vector{<:GACODE.FluxSolution}
end

"""
    ActorAnalyticalTurbulence(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates analytical turbulence models
"""
function ActorAnalyticalTurbulence(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorAnalyticalTurbulence(dd, act.ActorAnalyticalTurbulence; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorAnalyticalTurbulence(dd::IMAS.dd, par::FUSEparameters__ActorAnalyticalTurbulence; kw...)
    logging_actor_init(ActorAnalyticalTurbulence)
    par = OverrideParameters(par; kw...)
    return ActorAnalyticalTurbulence(dd, par, GACODE.FluxSolution[])
end

"""
    _step(actor::ActorAnalyticalTurbulence)

Runs analytical turbulent transport model on a vector of gridpoints
"""
function _step(actor::ActorAnalyticalTurbulence)
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]

    flux_solution = GACODE.FluxSolution(1.0, 1.0, 1.0, [1.0 for ion in cp1d.ion], 1.0)
    actor.flux_solutions = [flux_solution for irho in par.rho_transport]

    return actor
end

"""
    _finalize(actor::ActorAnalyticalTurbulence)

Writes results to dd.core_transport
"""
function _finalize(actor::ActorAnalyticalTurbulence)
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
