import NeoclassicalTransport
import GACODE

#= ================= =#
#  ActorNeoclassical  #
#= ================= =#
@actor_parameters_struct ActorNeoclassical{T} begin
    model::Switch{Symbol} = Switch{Symbol}([:changhinton, :neo, :hirshmansigmar], "-", "Neoclassical model to run"; default=:hirshmansigmar)
    neo_backend::Switch{Symbol} = Switch{Symbol}([:julia, :fortran], "-",
        "NEO backend for model=:neo: the native Julia port (run_neo_native, in process, threaded over the grid points, forward-AD capable), or the external Fortran NEO binary (run_neo)";
        default=:julia)
    collision_model::Entry{Int} = Entry{Int}("-",
        "NEO collision model: 1 Connor, 2 reduced Hirshman-Sigmar, 3 full Hirshman-Sigmar, 4 full linearized Fokker-Planck, 5 Fokker-Planck with ad-hoc field-particle terms";
        default=4)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.25:0.1:0.85)
end

mutable struct ActorNeoclassical{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorNeoclassical{P}}
    input_neos::Vector{<:NeoclassicalTransport.InputNEO}
    flux_solutions::Vector{GACODE.FluxSolution{D}}
    equilibrium_geometry::Union{NeoclassicalTransport.EquilibriumGeometry,Missing}
    neo_solutions::Vector{NeoclassicalTransport.NEOSolution{D}}
    neo_caches::Vector{NeoclassicalTransport.NEOFactorCache}
end

"""
    ActorNeoclassical(dd::IMAS.DD, act::ParametersAllActors; kw...)

Evaluates neoclassical (collisional) transport fluxes using established theoretical models.

Supported neoclassical models:
- `:changhinton`: Chang-Hinton model for ion heat transport in the banana/plateau regime
- `:neo`: Full drift-kinetic NEO code for comprehensive neoclassical transport including 
  bootstrap current, providing electron/ion energy, particle, and momentum fluxes.
  `neo_backend=:julia` (default) runs NeoclassicalTransport's native port in process
  (threaded over the grid points, usable with `jacobian_method=:forward_ad`) and keeps the
  full per-species `NEOSolution`s (bootstrap current, parallel flows, poloidal/toroidal
  velocities) in `actor.neo_solutions`; `neo_backend=:fortran` shells out to the GACODE
  NEO binary. `collision_model` selects the NEO collision operator for both backends
  (default 4, full Fokker-Planck).
- `:hirshmansigmar`: Hirshman-Sigmar analytical model for comprehensive neoclassical transport
  in various collisionality regimes

The models account for collisional effects, magnetic geometry, and trapped particle physics
to calculate transport coefficients. Results are stored in `dd.core_transport.model[:neoclassical]`
and include contributions to bootstrap current, thermal transport, and particle transport.
"""
function ActorNeoclassical(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorNeoclassical(dd, act.ActorNeoclassical; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorNeoclassical(dd::IMAS.DD{D}, par::FUSEparameters__ActorNeoclassical{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorNeoclassical)
    par = OverrideParameters(par; kw...)
    return ActorNeoclassical(dd, par, NeoclassicalTransport.InputNEO[], GACODE.FluxSolution{D}[], missing, NeoclassicalTransport.NEOSolution{D}[],
        NeoclassicalTransport.NEOFactorCache[])
end

"""
    neo_inputs(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, gridpoint_cps::AbstractVector{Int}, par)

The `InputNEO`s of `model=:neo` at the given `cp1d` grid points, with the actor's
`collision_model` applied. Generic in the `dd` number type, so the forward-AD
flux-matcher path can build Dual-valued inputs with it.
"""
function neo_inputs(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, gridpoint_cps::AbstractVector{Int}, par)
    input_neos = [NeoclassicalTransport.InputNEO(eqt, cp1d, i) for i in gridpoint_cps]
    for input_neo in input_neos
        input_neo.COLLISION_MODEL = par.collision_model
    end
    return input_neos
end

"""
    _step(actor::ActorNeoclassical)

Runs the selected neoclassical transport model to evaluate collisional fluxes on radial grid points.

For Chang-Hinton: Calculates ion heat transport using local parameters.
For NEO: Creates InputNEO structures and runs NEO for each grid point (the native Julia
solver as one threaded batch, or the Fortran binary through `asyncmap`).
For Hirshman-Sigmar: Uses cached equilibrium geometry and evaluates the analytical model
with local plasma parameters, providing comprehensive neoclassical transport coefficients.
"""
function _step(actor::ActorNeoclassical)
    par = actor.par
    dd = actor.dd

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    rho_cp = cp1d.grid.rho_tor_norm

    if par.model == :changhinton
        actor.flux_solutions = [NeoclassicalTransport.changhinton(eqt, cp1d, rho, 1) for rho in par.rho_transport]

    elseif par.model == :neo
        gridpoint_cps = [argmin_abs(rho_cp, rho) for rho in par.rho_transport]
        actor.input_neos = neo_inputs(eqt, cp1d, gridpoint_cps, par)
        if par.neo_backend == :julia
            # one factorization cache per grid point, kept across calls: flux-matcher
            # evaluations at nearby profiles reuse the previous factorization as a GMRES
            # preconditioner (refine=true), forward-AD passes at the same point reuse it
            # outright, and a new point on the same grid reuses the symbolic analysis
            if length(actor.neo_caches) != length(actor.input_neos)
                actor.neo_caches = [NeoclassicalTransport.NEOFactorCache(; refine=true) for _ in actor.input_neos]
            end
            actor.neo_solutions = NeoclassicalTransport.run_neo_native(actor.input_neos; caches=actor.neo_caches)
            actor.flux_solutions = GACODE.FluxSolution.(actor.neo_solutions)
        else
            actor.flux_solutions = asyncmap(input_neo -> NeoclassicalTransport.run_neo(input_neo), actor.input_neos)
        end

    elseif par.model == :hirshmansigmar
        gridpoint_cps = [argmin_abs(rho_cp, rho) for rho in par.rho_transport]
        if ismissing(actor.equilibrium_geometry) || actor.equilibrium_geometry.time != eqt.time
            actor.equilibrium_geometry = NeoclassicalTransport.get_equilibrium_geometry(eqt, cp1d)#, gridpoint_cps)
        end
        parameter_matrices = NeoclassicalTransport.get_plasma_profiles(eqt, cp1d)
        rho_s = GACODE.rho_s(cp1d, eqt)
        rmin = GACODE.r_min_core_profiles(eqt.profiles_1d, cp1d.grid.rho_tor_norm)
        actor.flux_solutions = map(gridpoint_cp -> NeoclassicalTransport.hirshmansigmar(gridpoint_cp, eqt, cp1d, parameter_matrices, actor.equilibrium_geometry; rho_s, rmin), gridpoint_cps)
    end

    return actor
end

"""
    _finalize(actor::ActorNeoclassical)

Writes neoclassical transport fluxes to `dd.core_transport.model[:neoclassical]`.

Sets the model identifier based on the selected model (Chang-Hinton, NEO, or Hirshman-Sigmar)
and converts flux results from GACODE format to IMAS format. The number and type of fluxes 
written depend on the model: Chang-Hinton provides only ion energy flux, while NEO and 
Hirshman-Sigmar provide comprehensive electron/ion energy, particle, and momentum fluxes.
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
        GACODE.flux_gacode_to_imas((:ion_energy_flux,), actor.flux_solutions, m1d, eqt, cp1d)

    elseif par.model == :neo
        model.identifier.name = par.neo_backend == :julia ? "NEO (Julia)" : "NEO"
        GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    elseif par.model == :hirshmansigmar
        model.identifier.name = "Hirshman-Sigmar"
        GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux), actor.flux_solutions, m1d, eqt, cp1d)
    end

    return actor
end
