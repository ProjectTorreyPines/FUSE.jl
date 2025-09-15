import TurbulentTransport
import TurbulentTransport: InputQLGYRO, InputCGYRO
import GACODE

#= =========== =#
#  ActorQLGYRO  #
#= =========== =#
@actor_parameters_struct ActorQLGYRO{T} begin
    ky::Entry{Float64} = Entry{Float64}("-", "Max ky"; default=1.6)
    nky::Entry{Int} = Entry{Int}("-", "Number of ky modes"; default=16)
    cpu_per_ky::Entry{Int} = Entry{Int}("-", "Number of cpus per ky"; default=1)
    kygrid_model::Entry{Int} = Entry{Int}("-", "TGLF ky grid model"; default=0)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat1, :sat2, :sat3], "-", "Saturation rule"; default=:sat1)
    n_field::Entry{Int} = Entry{Int}("-", "1:phi, 2:phi+apar, 3:phi+apar+bpar"; default=1, check=x -> @assert 1 <= x <= 3 "n_fields must be either 1,2 or 3")
    delta_t::Entry{Float64} = Entry{Float64}("-", "CGYRO step size "; default=0.005)
    max_time::Entry{Float64} = Entry{Float64}("-", "Max simulation time (a/cs)"; default=100.0)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute QLGYRO fluxes on"; default=0.25:0.1:0.85)
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
end

mutable struct ActorQLGYRO{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorQLGYRO{P}}
    input_qlgyros::Vector{InputQLGYRO}
    input_cgyros::Vector{InputCGYRO}
    flux_solutions::Vector{GACODE.FluxSolution{D}}
end

"""
    ActorQLGYRO(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates turbulent transport using the QLGYRO quasi-linear gyrokinetic model.

QLGYRO combines TGLF quasi-linear theory with CGYRO nonlinear gyrokinetic simulations 
to provide high-fidelity turbulent transport predictions. The model:

- Uses TGLF to identify the most unstable modes and calculate quasi-linear fluxes
- Runs CGYRO nonlinear gyrokinetic simulations to compute saturation levels
- Applies saturation rules to calibrate the quasi-linear transport predictions

Key parameters include ky spectral resolution, number of field components (electrostatic
vs electromagnetic), simulation time parameters, and saturation rules. The model provides
comprehensive electron/ion energy, particle, and momentum fluxes with gyrokinetic fidelity
while being more computationally efficient than full nonlinear simulations.
"""
function ActorQLGYRO(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorQLGYRO(dd, act.ActorQLGYRO; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorQLGYRO(dd::IMAS.dd{D}, par::FUSEparameters__ActorQLGYRO{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorQLGYRO)
    par = OverrideParameters(par; kw...)
    input_QLGYROs = Vector{InputQLGYRO}(undef, length(par.rho_transport))
    input_CGYROs = Vector{InputCGYRO}(undef, length(par.rho_transport))
    return ActorQLGYRO(dd, par, input_QLGYROs, input_CGYROs, GACODE.FluxSolution{D}[])
end

"""
    _step(actor::ActorQLGYRO)

Runs QLGYRO to evaluate gyrokinetic turbulent transport on radial grid points.

Creates InputQLGYRO and InputCGYRO input structures for each transport grid point,
configuring ky spectral parameters, field components, time stepping, and saturation rules.
Calls `TurbulentTransport.run_qlgyro` to execute the quasi-linear gyrokinetic calculation
combining TGLF linear analysis with CGYRO nonlinear saturation physics.
"""
function _step(actor::ActorQLGYRO)
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    ix_cp = [argmin_abs(cp1d.grid.rho_tor_norm, rho) for rho in par.rho_transport]

    for (k, gridpoint_cp) in enumerate(ix_cp)
        actor.input_qlgyros[k] = TurbulentTransport.InputQLGYRO()

        actor.input_qlgyros[k].N_PARALLEL = par.nky * par.cpu_per_ky
        actor.input_qlgyros[k].N_RUNS = 1
        actor.input_qlgyros[k].GAMMA_E = 0.0
        actor.input_qlgyros[k].CODE = -1
        actor.input_qlgyros[k].NKY = par.nky
        actor.input_qlgyros[k].KYGRID_MODEL = par.kygrid_model
        actor.input_qlgyros[k].KY = par.ky

        if par.sat_rule == :sat2
            actor.input_qlgyros[k].SAT_RULE = 2
        elseif par.sat_rule == :sat1
            actor.input_qlgyros[k].SAT_RULE = 1
        end

        actor.input_cgyros[k] = InputCGYRO(dd, gridpoint_cp, par.lump_ions)
        actor.input_cgyros[k].N_FIELD = par.n_field
        actor.input_cgyros[k].DELTA_T = par.delta_t
        actor.input_cgyros[k].MAX_TIME = par.max_time

        actor.input_cgyros[k].DELTA_T_METHOD = 1
        actor.input_cgyros[k].FREQ_TOL = 0.01
    end

    actor.flux_solutions = TurbulentTransport.run_qlgyro(actor.input_qlgyros, actor.input_cgyros)

    return actor
end

"""
    _finalize(actor::ActorQLGYRO)

Writes QLGYRO transport fluxes to `dd.core_transport.model[:anomalous]`.

Sets the model identifier to "QLGYRO" and converts flux results from GACODE format 
to IMAS format for electron/ion energy, particle, and momentum fluxes using
`GACODE.flux_gacode_to_imas`.
"""
function _finalize(actor::ActorQLGYRO)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = "QLGYRO"
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end

function Base.show(io::IO, ::MIME"text/plain", input::InputQLGYRO)
    for field_name in fieldnames(typeof(input))
        println(io, " $field_name = $(getfield(input,field_name))")
    end
end
