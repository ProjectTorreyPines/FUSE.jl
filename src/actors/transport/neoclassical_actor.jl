import TGLFNN
import TAUENN
import NEO

#= ===================== =#
#  ActorNeoclassical      #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorNeoclassical{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    neoclassical_model::Switch{Symbol} = Switch{Symbol}([:changhinton, :neo], "-", "Neoclassical model to run"; default=:changhinton)
    rho_transport::Entry{AbstractVector{<:T}} = Entry{AbstractVector{<:T}}("-", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.2:0.1:0.8)
end

mutable struct ActorNeoclassical{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorNeoclassical{P}
    flux_solutions::Vector{<:TGLFNN.flux_solution}
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
    _step(actor::ActorNeoclassical)

Runs ActorNeoclassical to evaluate the neoclassical transport flux on a vector of gridpoints
"""
function _step(actor::ActorNeoclassical)
    par = actor.par
    dd = actor.dd
    model = resize!(dd.core_transport.model, :neoclassical; wipe=false)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport
    cp1d = dd.core_profiles.profiles_1d[]

    if par.neoclassical_model == :changhinton
        model.identifier.name = "Chang-Hinton"
        eqt = dd.equilibrium.time_slice[]
        actor.flux_solutions = [TAUENN.neoclassical_changhinton(eqt, cp1d, rho, 1) for rho in par.rho_transport]
    elseif par.neoclassical_model == :neo 
        model.identifier.name = "NEO"
        rho_cp = cp1d.grid.rho_tor_norm
        gridpoint_cp = [argmin(abs.(rho_cp .- rho)) for rho in par.rho_transport]

        n_species = length(cp1d.ion)+1
        
        for i in gridpoint_cp
            neo_solution = NEO.run_neo(NEO.InputNEO(dd, i))
           
            energy_flux_electrons = getfield(neo_solution, Symbol("ENERGY_FLUX_$(n_species)")) # electrons are always the last species in the list
            total_ion_energy_flux = -energy_flux_electrons # exclude electron energy flux from total ion energy flux
            particle_flux_electrons = getfield(neo_solution, Symbol("PARTICLE_FLUX_$(n_species)"))

            for field in fieldnames(NEO.flux_solution)
                if ismissing(getfield(neo_solution, field))
                    setfield!(neo_solution, field, 0.0)
                end

                if occursin("ENERGY", string(field))
                    total_ion_energy_flux += getfield(neo_solution, field) 
                end
            end

            solution = TGLFNN.flux_solution(particle_flux_electrons, 0.0, energy_flux_electrons, total_ion_energy_flux)
            push!(actor.flux_solutions, solution)
        end
  
    else
        error("$(par.neoclassical_model) is not implemented")
    end

    return actor
end

"""
    _finalize(actor::ActorNeoclassical)

Writes ActorNeoclassical results to dd.core_transport
"""
function _finalize(actor::ActorNeoclassical)
    par = actor.par
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = findfirst(:neoclassical, actor.dd.core_transport.model)
    m1d = model.profiles_1d[]
    m1d.total_ion_energy.flux = zeros(length(par.rho_transport))
    m1d.electrons.particles.flux = zeros(length(par.rho_transport))
    m1d.electrons.energy.flux = zeros(length(par.rho_transport))

    for (neoclassical_idx, rho) in enumerate(par.rho_transport)
        rho_transp_idx = findfirst(i -> i == rho, m1d.grid_flux.rho_tor_norm)
        rho_cp_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho))
        m1d.total_ion_energy.flux[rho_transp_idx] = actor.flux_solutions[neoclassical_idx].ENERGY_FLUX_i * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2

        if par.neoclassical_model == :neo
            m1d.electrons.particles.flux[rho_transp_idx] = actor.flux_solutions[neoclassical_idx].PARTICLE_FLUX_e * IMAS.gyrobohm_particle_flux(cp1d, eqt)[rho_cp_idx]
            m1d.electrons.energy.flux[rho_transp_idx] = actor.flux_solutions[neoclassical_idx].ENERGY_FLUX_e * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx]
        end
    end

    return actor
end
