import TGLFNN

#= ============= =#
#  ActorTGLF      #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorTGLF{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    nn::Entry{Bool} = Entry{Bool}("-", "Use TGLF-NN"; default=true)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat0quench, :sat1, :sat1geo, :sat2], "-", "Saturation rule"; default=:sat1)
    electromagnetic::Entry{Bool} = Entry{Bool}("-", "Electromagnetic or electrostatic"; default=true)
    rho_transport::Entry{AbstractVector{<:T}} = Entry{AbstractVector{<:T}}("-", "rho_tor_norm values to compute tglf fluxes on"; default=0.2:0.1:0.8)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
end

mutable struct ActorTGLF{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTGLF{P}
    input_tglfs::Vector{<:TGLFNN.InputTGLF}
    flux_solutions::Vector{<:TGLFNN.flux_solution}
end

"""
    ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the TGLF predicted turbulence
"""
function ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTGLF(dd, act.ActorTGLF; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTGLF(dd::IMAS.dd, par::FUSEparameters__ActorTGLF; kw...)
    par = par(kw...)
    input_tglfs = Vector{TGLFNN.InputTGLF}(undef, length(par.rho_transport))
    return ActorTGLF(dd, par, input_tglfs, TGLFNN.flux_solution[])
end

"""
    step(actor::ActorTGLF)

Runs TGLF actor to evaluate the turbulence flux on a vector of gridpoints
"""
function _step(actor::ActorTGLF)
    par = actor.par
    dd = actor.dd

    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    ix_eq = [argmin(abs.(eq1d.rho_tor_norm .- rho)) for rho in par.rho_transport]
    ix_cp = [argmin(abs.(cp1d.grid.rho_tor_norm .- rho)) for rho in par.rho_transport]
    for (k, (gridpoint_eq, gridpoint_cp)) in enumerate(zip(ix_eq, ix_cp))
        actor.input_tglfs[k] = TGLFNN.InputTGLF(dd, gridpoint_eq, gridpoint_cp, par.sat_rule, par.electromagnetic)
    end

    model_filename = string(par.sat_rule) * "_" * (par.electromagnetic ? "em" : "es")

    model = resize!(dd.core_transport.model, :anomalous)
    model.identifier.name = (par.nn ? "TGLF-NN" : "TGLF") * " " * model_filename
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    if par.nn
        actor.flux_solutions = map(input_tglf -> TGLFNN.run_tglfnn(input_tglf; par.warn_nn_train_bounds, model_filename), actor.input_tglfs)
    else
        actor.flux_solutions = asyncmap(input_tglf -> TGLFNN.run_tglf(input_tglf), actor.input_tglfs)
    end

    return actor
end

"""
    finalize(actor::ActorTGLF)

Writes results to dd.core_transport
"""
function _finalize(actor::ActorTGLF)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    model = findfirst(:anomalous, actor.dd.core_transport.model)
    m1d = model.profiles_1d[]
    m1d.electrons.energy.flux = zeros(length(par.rho_transport))
    m1d.total_ion_energy.flux = zeros(length(par.rho_transport))
    m1d.electrons.particles.flux = zeros(length(par.rho_transport))
    m1d.momentum_tor.flux = zeros(length(par.rho_transport))
    for (tglf_idx, rho) in enumerate(par.rho_transport)
        rho_transp_idx = findfirst(i -> i == rho, m1d.grid_flux.rho_tor_norm)
        rho_cp_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho))
        m1d.electrons.energy.flux[rho_transp_idx] = actor.flux_solutions[tglf_idx].ENERGY_FLUX_e * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2
        m1d.total_ion_energy.flux[rho_transp_idx] = actor.flux_solutions[tglf_idx].ENERGY_FLUX_i * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2
        m1d.electrons.particles.flux[rho_transp_idx] = actor.flux_solutions[tglf_idx].PARTICLE_FLUX_e * IMAS.gyrobohm_particle_flux(cp1d, eqt)[rho_cp_idx] # 1 / m^2 / s
        m1d.momentum_tor.flux[rho_transp_idx] = actor.flux_solutions[tglf_idx].STRESS_TOR_i * IMAS.gyrobohm_momentum_flux(cp1d, eqt)[rho_cp_idx] #
    end
    return actor
end
