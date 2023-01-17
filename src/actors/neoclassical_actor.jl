import TGLFNN
import TAUENN

#= ===================== =#
#  ActorNeoclassical      #
#= ===================== =#
mutable struct ActorNeoclassical <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    flux_solutions::AbstractVector{<:TGLFNN.flux_solution}
end

Base.@kwdef struct FUSEparameters__ActorNeoclassical{T} <: ParametersActor where {T<:Real}
    neoclassical_model = Switch(Symbol, [:changhinton], "", "Neoclassical model to run"; default=:changhinton)
    rho_transport = Entry(AbstractVector{<:Real}, "", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.2:0.1:0.8)
end

"""
    ActorNeoclassical(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the neoclassical transport fluxes
"""
function ActorNeoclassical(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorNeoclassical(kw...)
    actor = ActorNeoclassical(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorNeoclassical(dd::IMAS.dd, par::ParametersActor; kw...)
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
    neoclassical_index = IMAS.name_2_index(dd.core_transport.model)[:neoclassical]
    model = resize!(dd.core_transport.model, "identifier.index" => neoclassical_index)
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
function finalize(actor::ActorNeoclassical)
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
