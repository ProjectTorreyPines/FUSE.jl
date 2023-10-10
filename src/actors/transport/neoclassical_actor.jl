import NEO
#= ===================== =#
#  ActorNeoclassical      #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorNeoclassical{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch{Symbol}([:changhinton, :neo], "-", "Neoclassical model to run"; default=:changhinton)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.2:0.05:0.85)
end

mutable struct ActorNeoclassical{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorNeoclassical{P}
    input_neos::Union{Vector{<:NEO.InputNEO},Missing}
    flux_solutions::Vector{<:IMAS.flux_solution}
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
    if par.model == :neo
        input_neos = Vector{NEO.InputNEO}(undef, length(par.rho_transport))
    else
        input_neos = missing
    end
    return ActorNeoclassical(dd, par, input_neos, IMAS.flux_solution[])
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

    if par.model == :changhinton
        model.identifier.name = "Chang-Hinton"
        eqt = dd.equilibrium.time_slice[]
        actor.flux_solutions = [NEO.changhinton(eqt, cp1d, rho, 1) for rho in par.rho_transport]
    elseif par.model == :neo
        model.identifier.name = "NEO"
        rho_cp = cp1d.grid.rho_tor_norm
        gridpoint_cp = [argmin(abs.(rho_cp .- rho)) for rho in par.rho_transport]
        for (idx, i) in enumerate(gridpoint_cp)
            actor.input_neos[idx] = NEO.InputNEO(dd, i)
        end
        actor.flux_solutions = asyncmap(input_neo -> NEO.run_neo(input_neo), actor.input_neos)
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

    if par.model == :changhinton
        IMAS.flux_gacode_to_fuse([:ion_energy_flux], actor.flux_solutions, m1d, eqt, cp1d)
    elseif par.model == :neo
        IMAS.flux_gacode_to_fuse([:ion_energy_flux, :electron_energy_flux, :electron_particle_flux], actor.flux_solutions, m1d, eqt, cp1d)
    end

    return actor
end
