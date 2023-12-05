import NEO

#= ================= =#
#  ActorNeoclassical  #
#= ================= =#
Base.@kwdef mutable struct FUSEparameters__ActorNeoclassical{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch{Symbol}([:changhinton, :neo, :hirshmansigmar], "-", "Neoclassical model to run"; default=:hirshmansigmar)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.25:0.1:0.85)
end

mutable struct ActorNeoclassical{D,P} <: PlasmaAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorNeoclassical{P}
    input_neos::Vector{<:NEO.InputNEO}
    flux_solutions::Vector{<:IMAS.flux_solution}
    equilibrium_geometry::Union{NEO.equilibrium_geometry,Missing}
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
    return ActorNeoclassical(dd, par, NEO.InputNEO[], IMAS.flux_solution[], missing)
end

"""
    _step(actor::ActorNeoclassical)

Runs ActorNeoclassical to evaluate the neoclassical transport flux on a vector of gridpoints
"""
function _step(actor::ActorNeoclassical)
    par = actor.par
    dd = actor.dd

    cp1d = dd.core_profiles.profiles_1d[]
    rho_cp = cp1d.grid.rho_tor_norm

    if par.model == :changhinton
        eqt = dd.equilibrium.time_slice[]
        actor.flux_solutions = [NEO.changhinton(eqt, cp1d, rho, 1) for rho in par.rho_transport]

    elseif par.model == :neo
        gridpoint_cps = [argmin(abs.(rho_cp .- rho)) for rho in par.rho_transport]
        actor.input_neos = [NEO.InputNEO(dd, i) for (idx, i) in enumerate(gridpoint_cps)]
        actor.flux_solutions = asyncmap(input_neo -> NEO.run_neo(input_neo), actor.input_neos)

    elseif par.model == :hirshmansigmar
        gridpoint_cps = [argmin(abs.(rho_cp .- rho)) for rho in par.rho_transport]
        if ismissing(actor.equilibrium_geometry)
            actor.equilibrium_geometry = NEO.get_equilibrium_parameters(actor.dd)
        end
        parameter_matrices = NEO.get_ion_electron_parameters(dd)
        actor.flux_solutions = map(gridpoint_cp -> NEO.hirshmansigmar(gridpoint_cp, dd, parameter_matrices, actor.equilibrium_geometry), gridpoint_cps)
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

    model = resize!(dd.core_transport.model, :neoclassical; wipe=false)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    if par.model == :changhinton
        model.identifier.name = "Chang-Hinton"
        IMAS.flux_gacode_to_fuse([:ion_energy_flux], actor.flux_solutions, m1d, eqt, cp1d)

    elseif par.model == :neo
        model.identifier.name = "NEO"
        IMAS.flux_gacode_to_fuse([:ion_energy_flux, :electron_energy_flux, :electron_particle_flux], actor.flux_solutions, m1d, eqt, cp1d)

    elseif par.model == :hirshmansigmar
        model.identifier.name = "Hirshman-Sigmar"
        IMAS.flux_gacode_to_fuse([:ion_energy_flux, :electron_energy_flux, :electron_particle_flux], actor.flux_solutions, m1d, eqt, cp1d)
    end

    return actor
end
