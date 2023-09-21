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
    user_specified_model::Entry{String} = Entry{String}("-", "Use a user specified TGLF-NN model stored in TGLFNN/models"; default = "")
    rho_transport::Entry{AbstractVector{<:T}} = Entry{AbstractVector{<:T}}("-", "rho_tor_norm values to compute tglf fluxes on"; default=0.2:0.1:0.8)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
end

mutable struct ActorTGLF{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTGLF{P}
    input_tglfs::Vector{<:TGLFNN.InputTGLF}
    flux_solutions::Vector{<:IMAS.flux_solution}
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
    return ActorTGLF(dd, par, input_tglfs, IMAS.flux_solution[])
end

"""
    _step(actor::ActorTGLF)

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

    if !isempty(par.user_specified_model)
        model_filename = par.user_specified_model
    else
        model_filename = string(par.sat_rule) * "_" * (par.electromagnetic ? "em" : "es")
        model_filename *= "_d3d" # will be changed to FPP soon
    end


    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = (par.nn ? "TGLF-NN" : "TGLF") * " " * model_filename
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    if par.nn
        actor.flux_solutions = TGLFNN.run_tglfnn(actor.input_tglfs; par.warn_nn_train_bounds, model_filename)
    else
        actor.flux_solutions = TGLFNN.run_tglf(actor.input_tglfs)
    end

    return actor
end

"""
    _finalize(actor::ActorTGLF)

Writes results to dd.core_transport
"""
function _finalize(actor::ActorTGLF)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    model = findfirst(:anomalous, actor.dd.core_transport.model)
    m1d = model.profiles_1d[]

    IMAS.flux_gacode_to_fuse([:ion_energy_flux, :electron_energy_flux, :electron_particle_flux, :momentum_flux], actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end
