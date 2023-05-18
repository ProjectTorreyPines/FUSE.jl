import TGLFNN
import TAUENN

#= ===================== =#
#  ActorNeoclassical      #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorNeoclassical{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    neoclassical_model::Switch{Symbol} = Switch{Symbol}([:changhinton], "-", "Neoclassical model to run"; default=:changhinton)
    rho_transport::Entry{AbstractVector{<:T}} = Entry{AbstractVector{<:T}}("-", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.2:0.1:0.8)
end

mutable struct ActorNeoclassical <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorNeoclassical
    flux_solutions::AbstractVector{<:TGLFNN.flux_solution}
end

"""
    ActorNeoclassical(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the neoclassical transport fluxes
"""
function ActorNeoclassical(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorNeoclassical(dd, act.ActorNeoclassical; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorNeoclassical(dd::IMAS.dd, par::FUSEparameters__ActorNeoclassical; kw...)
    par = par(kw...)
    return ActorNeoclassical(dd, par, TGLFNN.flux_solution[])
end

"""
    step(actor::ActorNeoclassical)

Runs ActorNeoclassical to evaluate the neoclassical transport flux on a vector of gridpoints
"""
function _step(actor::ActorNeoclassical)
    par = actor.par
    dd = actor.dd
    model = resize!(dd.core_transport.model, :neoclassical)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    if par.neoclassical_model == :changhinton
        model.identifier.name = "Chang-Hinton"
        eqt = dd.equilibrium.time_slice[]
        cp1d = dd.core_profiles.profiles_1d[]
        actor.flux_solutions = [TAUENN.neoclassical_changhinton(eqt, cp1d, rho, 1) for rho in par.rho_transport]
    else
        error("$(par.neoclassical_model) is not implemented")
    end

    return actor
end

"""
    finalize(actor::ActorNeoclassical)

Writes ActorNeoclassical results to dd.core_transport
"""
function _finalize(actor::ActorNeoclassical)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = findfirst(:neoclassical, actor.dd.core_transport.model)
    m1d = model.profiles_1d[]
    m1d.total_ion_energy.flux = zeros(length(actor.par.rho_transport))
    for (neoclassical_idx, rho) in enumerate(actor.par.rho_transport)
        rho_transp_idx = findfirst(i -> i == rho, m1d.grid_flux.rho_tor_norm)
        rho_cp_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho))
        m1d.total_ion_energy.flux[rho_transp_idx] = actor.flux_solutions[neoclassical_idx].ENERGY_FLUX_i * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2
    end

    return actor
end
