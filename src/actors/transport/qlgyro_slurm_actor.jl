import TurbulentTransport
import TurbulentTransport: InputQLGYRO, InputCGYRO
import GACODE

#= ================ =#
#  ActorQLGYROSlurm  #
#= ================ =#
@actor_parameters_struct ActorQLGYROSlurm{T} begin
    ky::Entry{Float64} = Entry{Float64}("-", "Max ky"; default=1.6)
    nky::Entry{Int} = Entry{Int}("-", "Number of ky modes"; default=16)
    kygrid_model::Entry{Int} = Entry{Int}("-", "TGLF ky grid model"; default=0)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat1, :sat2, :sat3], "-", "Saturation rule"; default=:sat1)
    n_field::Entry{Int} = Entry{Int}("-", "1:phi, 2:phi+apar, 3:phi+apar+bpar"; default=1, check=x -> @assert 1 <= x <= 3 "n_fields must be either 1,2 or 3")
    delta_t::Entry{Float64} = Entry{Float64}("-", "CGYRO step size "; default=0.005)
    max_time::Entry{Float64} = Entry{Float64}("-", "Max simulation time (a/cs)"; default=100.0)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute QLGYRO fluxes on"; default=0.25:0.1:0.85)
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
    MXH_modes::Entry{Int} = Entry{Int}("-", "Number of MXH harmonic modes (1 = standard Miller geometry)"; default=1)
    basedir::Entry{String} = Entry{String}("-", "Root directory for CGYRO run trees (one subdir per rho); defaults to the shared m3739 results area, falling back to \$PSCRATCH"; default="")
    gpu::Entry{Bool} = Entry{Bool}("-", "Run CGYRO on the GPU queue"; default=true)
    n_mpi::Entry{Int} = Entry{Int}("-", "Number of MPI ranks per CGYRO ky run"; default=32)
    n_omp::Entry{Int} = Entry{Int}("-", "Number of OpenMP threads per CGYRO ky run"; default=4)
    walltime::Entry{String} = Entry{String}("-", "SLURM walltime per CGYRO ky job"; default="00:15:00")
    repo::Entry{String} = Entry{String}("-", "NERSC project allocation to submit CGYRO jobs under"; default="m3739_g")
    qos::Switch{Symbol} = Switch{Symbol}([:regular, :premium, :debug], "-", "SLURM QOS for CGYRO jobs"; default=:regular)
    poll_interval::Entry{Int} = Entry{Int}("-", "Seconds between convergence checks while waiting for CGYRO jobs"; default=60)
end

mutable struct ActorQLGYROSlurm{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorQLGYROSlurm{P}}
    input_qlgyros::Vector{InputQLGYRO}
    input_cgyros::Vector{InputCGYRO}
    flux_solutions::Vector{GACODE.FluxSolution{D}}
end

"""
    ActorQLGYROSlurm(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates turbulent transport using the QLGYRO quasi-linear gyrokinetic model, submitting
CGYRO's linear ky runs to SLURM under a specific NERSC allocation rather than running the
monolithic `qlgyro` executable in-process (as `ActorQLGYRO` does).

QLGYRO combines TGLF quasi-linear theory with CGYRO nonlinear gyrokinetic simulations
to provide high-fidelity turbulent transport predictions. The model:

- Uses TGLF to identify the most unstable modes and calculate quasi-linear fluxes
- Runs CGYRO nonlinear gyrokinetic simulations to compute saturation levels
- Applies saturation rules to calibrate the quasi-linear transport predictions

Each radial point's `(InputCGYRO, InputQLGYRO)` pair drives its own linear CGYRO ky
scan, submitted as SLURM batch jobs (via `repo`/`qos`/`walltime`/`gpu`/`n_mpi`/`n_omp`).
This lets CGYRO run under a NERSC allocation distinct from whatever resources FUSE
itself is running on, and is what's needed when this actor is used as the
`turbulence_model` inside `ActorFluxMatcher`/`ActorFluxCalculator`, since the flux
evaluation still returns synchronously (blocks until all ky's for a radius have
converged) so the nonlinear solve driving `ActorFluxMatcher` can proceed.
"""
function ActorQLGYROSlurm(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorQLGYROSlurm(dd, act.ActorQLGYROSlurm; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorQLGYROSlurm(dd::IMAS.dd{D}, par::FUSEparameters__ActorQLGYROSlurm{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorQLGYROSlurm)
    par = OverrideParameters(par; kw...)
    input_QLGYROs = Vector{InputQLGYRO}(undef, length(par.rho_transport))
    input_CGYROs = Vector{InputCGYRO}(undef, length(par.rho_transport))
    return ActorQLGYROSlurm(dd, par, input_QLGYROs, input_CGYROs, GACODE.FluxSolution{D}[])
end

"""
    _step(actor::ActorQLGYROSlurm)

Runs QLGYRO to evaluate gyrokinetic turbulent transport on radial grid points.

Creates `InputQLGYRO`/`InputCGYRO` input structures for each transport grid point
(same construction `ActorQLGYRO`'s legacy in-process path uses), then for each radius
calls the modular `TurbulentTransport.run_qlgyro(::InputCGYRO, ::InputQLGYRO; ...)`
path, which submits a linear CGYRO run per ky as a SLURM batch job under `repo`/`qos`
and blocks until all ky's for that radius have converged, applying TJLF saturation
rules to return a `GACODE.FluxSolution`.
"""
function _step(actor::ActorQLGYROSlurm{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    ix_cp = [argmin_abs(cp1d.grid.rho_tor_norm, rho) for rho in par.rho_transport]

    actor.flux_solutions = Vector{GACODE.FluxSolution{D}}(undef, length(par.rho_transport))

    for (k, gridpoint_cp) in enumerate(ix_cp)
        input_qlgyro = TurbulentTransport.InputQLGYRO()
        input_qlgyro.NKY = par.nky
        input_qlgyro.KYGRID_MODEL = par.kygrid_model
        input_qlgyro.KY = par.ky
        if par.sat_rule == :sat1
            input_qlgyro.SAT_RULE = 1
        elseif par.sat_rule == :sat2
            input_qlgyro.SAT_RULE = 2
        elseif par.sat_rule == :sat3
            input_qlgyro.SAT_RULE = 3
        end
        actor.input_qlgyros[k] = input_qlgyro

        input_cgyro = InputCGYRO(dd, gridpoint_cp, par.lump_ions; MXH_modes=par.MXH_modes)
        input_cgyro.N_FIELD = par.n_field
        input_cgyro.DELTA_T = par.delta_t
        input_cgyro.MAX_TIME = par.max_time
        input_cgyro.DELTA_T_METHOD = 1
        input_cgyro.FREQ_TOL = 0.01
        actor.input_cgyros[k] = input_cgyro

        rho = par.rho_transport[k]
        rundir = isempty(par.basedir) ? "" : joinpath(par.basedir, "rho_$(round(rho; digits=4))")

        sol = TurbulentTransport.run_qlgyro(input_cgyro, input_qlgyro;
            basedir=rundir,
            gpu=par.gpu,
            n_mpi=par.n_mpi,
            n_omp=par.n_omp,
            walltime=par.walltime,
            repo=par.repo,
            qos=String(par.qos),
            wait_for_completion=true,
            poll_interval=par.poll_interval)

        actor.flux_solutions[k] = GACODE.FluxSolution{D}(
            D(sol.ENERGY_FLUX_e),
            D(sol.ENERGY_FLUX_i),
            D(sol.PARTICLE_FLUX_e),
            D.(sol.PARTICLE_FLUX_i),
            D(sol.STRESS_TOR_i))
    end

    return actor
end

"""
    _finalize(actor::ActorQLGYROSlurm)

Writes QLGYRO transport fluxes to `dd.core_transport.model[:anomalous]`.

Sets the model identifier to "QLGYRO" and converts flux results from GACODE format
to IMAS format for electron/ion energy, particle, and momentum fluxes using
`GACODE.flux_gacode_to_imas`.
"""
function _finalize(actor::ActorQLGYROSlurm)
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
