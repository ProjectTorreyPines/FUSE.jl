import GACODE
import JSON

#= ====================== =#
#  ActorFluxMatcherGACODE  #
#= ====================== =#
@actor_parameters_struct ActorFluxMatcherGACODE{T} begin
    executable::Entry{String} = Entry{String}("-", "Command used to launch the external flux-matcher (may include arguments, split on whitespace)")
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "ρ transport grid the external code flux-matches on"; default=0.25:0.1:0.85)
    predicted_channels::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "Channels to flux-match (subset of `:te`, `:ti`, `:ne`)"; default=[:te, :ti, :ne])
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat1, :sat2, :sat3], "-", "TGLF saturation rule used by the external code"; default=:sat2)
    max_iterations::Entry{Int} = Entry{Int}("-", "Maximum number of optimization iterations for the external flux-matcher"; default=20)
    input_gacode_filename::Entry{String} = Entry{String}("-", "Name of the input.gacode file written for the external code"; default="input.gacode")
    output_gacode_filename::Entry{String} =
        Entry{String}("-", "Name of the input.gacode file the external code writes with flux-matched profiles"; default="input.gacode.new")
    save_dir::Entry{String} = Entry{String}("-", "Working directory for the run (empty creates a temporary directory)"; default="")
    clear_workdir::Entry{Bool} = Entry{Bool}("-", "Remove the working directory after the run (only when it was auto-created)"; default=true)
    save_fluxes::Entry{Bool} = Entry{Bool}("-", "Read the external code's turbulent/neoclassical fluxes and store them in `dd.core_transport`"; default=true)
    evolve_rotation::Switch{Symbol} = Switch{Symbol}(
        [:flux_match, :fixed],
        "-",
        "Rotation evolution: `:flux_match` sets rotation from the external code's `w0`, `:fixed` leaves the existing rotation profile alone";
        default=:fixed)
end

mutable struct ActorFluxMatcherGACODE{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorFluxMatcherGACODE{P}}
    input_gacode_out::Union{GACODE.InputGACODE,Missing}
    flux_data::Dict{Symbol,Any}  # parsed `fluxes_turb.json` / `fluxes_neoc.json` keyed by :turb / :neoc
end

"""
    ActorFluxMatcherGACODE(dd::IMAS.dd, act::ParametersAllActors; kw...)

Flux-matches core plasma profiles using an external code that communicates through `input.gacode` files.

The actor:

 1. Builds an `input.gacode` from the current `dd` (`GACODE.InputGACODE(dd)`) and writes it to the working directory.
 2. Runs the external executable in that directory; the code reads the `input.gacode` and is expected to write a
    new `input.gacode` file (`output_gacode_filename`) containing the flux-matched profiles.
 3. Loads the returned file and updates `dd.core_profiles` with the new profiles.

Updated channels:

  - Electron temperature (`te`)
  - Ion temperature (`ti`, applied to all thermal ions)
  - Electron density (`ne`)
  - Rotation (`w0`) when `evolve_rotation == :flux_match`

The external code is responsible for producing the full radial profile (core + edge); profiles are interpolated
from the returned grid onto the current `core_profiles` grid.
"""
function ActorFluxMatcherGACODE(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFluxMatcherGACODE(dd, act.ActorFluxMatcherGACODE; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFluxMatcherGACODE(dd::IMAS.dd{D}, par::FUSEparameters__ActorFluxMatcherGACODE{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorFluxMatcherGACODE)
    par = OverrideParameters(par; kw...)
    return ActorFluxMatcherGACODE(dd, par, missing, Dict{Symbol,Any}())
end

"""
    _step(actor::ActorFluxMatcherGACODE)

Writes the current `dd` as an `input.gacode`, runs the external flux-matcher, and loads the returned `input.gacode`.
"""
function _step(actor::ActorFluxMatcherGACODE)
    dd = actor.dd
    par = actor.par

    @assert !isempty(strip(par.executable)) "ActorFluxMatcherGACODE: `executable` must be set to the external flux-matcher command"

    # Working directory: reuse `save_dir` if given, otherwise create a temporary one
    created_workdir = isempty(par.save_dir)
    run_dir = created_workdir ? mktempdir() : abspath(par.save_dir)
    mkpath(run_dir)
    par.save_dir = run_dir
    @info "ActorFluxMatcherGACODE running in $run_dir"

    # Write the input.gacode for the external code
    input_gacode = GACODE.InputGACODE(dd)
    input_path = joinpath(run_dir, par.input_gacode_filename)
    GACODE.save(input_gacode, input_path)

    # Build the command: the user-supplied `executable` plus the flux-match settings
    # (rho grid, channels, saturation rule, iterations, file names) as CLI flags so
    # the run is fully controlled from the FUSE side.
    flags = [
        "--input", par.input_gacode_filename,
        "--output", par.output_gacode_filename,
        "--predicted-rho", join(string.(collect(par.rho_transport)), ","),
        "--predicted-channels", join(string.(par.predicted_channels), ","),
        "--code-settings", uppercase(string(par.sat_rule)),
        "--max-iterations", string(par.max_iterations)
    ]
    cmd = Cmd(vcat(String.(split(strip(par.executable))), flags))

    # Run the external code in the working directory
    old_dir = pwd()
    try
        cd(run_dir)
        Base.run(cmd)  # qualify: FUSE defines its own `run` for studies, which shadows Base.run
    finally
        cd(old_dir)
    end

    # Load the flux-matched profiles returned by the external code
    output_path = joinpath(run_dir, par.output_gacode_filename)
    @assert isfile(output_path) "ActorFluxMatcherGACODE: expected output file `$output_path` was not produced by `$(par.executable)`"
    actor.input_gacode_out = GACODE.load(output_path)

    # Read the turbulent/neoclassical fluxes (gyro-Bohm units) the external code wrote for
    # the best evaluation. Must happen before any cleanup of the working directory.
    empty!(actor.flux_data)
    if par.save_fluxes
        for (key, fname) in ((:turb, "fluxes_turb.json"), (:neoc, "fluxes_neoc.json"))
            fpath = joinpath(run_dir, fname)
            if isfile(fpath)
                actor.flux_data[key] = JSON.parsefile(fpath)
            else
                @warn "ActorFluxMatcherGACODE: `$fname` not found in $run_dir; skipping $key fluxes"
            end
        end
    end

    # Only auto-clean a directory we created ourselves
    if created_workdir && par.clear_workdir
        rm(run_dir; force=true, recursive=true)
    end

    return actor
end

"""
    _finalize(actor::ActorFluxMatcherGACODE)

Updates `dd.core_profiles` from the flux-matched `input.gacode` returned by the external code.

Profiles are interpolated from the returned `rho_tor_norm` grid onto the current `core_profiles` grid.
Temperatures are stored in eV (gacode is keV), densities in m⁻³ (gacode is 10¹⁹ m⁻³), and rotation `w0` maps
directly to `rotation_frequency_tor_sonic`.
"""
function _finalize(actor::ActorFluxMatcherGACODE)
    dd = actor.dd
    par = actor.par
    ig = actor.input_gacode_out

    @assert !ismissing(ig) "ActorFluxMatcherGACODE: no output loaded; run `step` before `finalize`"

    cp1d = dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm
    rho_out = ig.rho

    # Electron temperature (keV → eV)
    cp1d.electrons.temperature = IMAS.interp1d(rho_out, ig.te .* 1e3).(rho)

    # Electron density (10¹⁹ m⁻³ → m⁻³)
    cp1d.electrons.density_thermal = IMAS.interp1d(rho_out, ig.ne .* 1e19).(rho)
    IMAS.unfreeze!(cp1d.electrons, :density)

    # Ion temperature: use the first thermal ion row and apply to all ions (keV → eV)
    Ti_new = IMAS.interp1d(rho_out, ig.ti[1, :] .* 1e3).(rho)
    for ion in cp1d.ion
        ion.temperature = Ti_new
    end
    IMAS.unfreeze!(cp1d, :t_i_average)

    # Rotation
    if par.evolve_rotation == :flux_match && !ismissing(ig.w0)
        cp1d.rotation_frequency_tor_sonic = IMAS.interp1d(rho_out, ig.w0).(rho)
    end

    IMAS.intrinsic_sources!(dd)

    # Store the external code's fluxes in dd.core_transport (same GB→IMAS path as ActorQLGYRO)
    if par.save_fluxes && !isempty(actor.flux_data)
        write_gacode_fluxes!(actor)
    end

    return actor
end

"""
    write_gacode_fluxes!(actor::ActorFluxMatcherGACODE)

Convert the external code's gyro-Bohm fluxes (`fluxes_turb.json` / `fluxes_neoc.json`) into
IMAS physical units and store them in `dd.core_transport.model[:anomalous]` (turbulent) and
`[:neoclassical]`, reusing `GACODE.flux_gacode_to_imas` (the same conversion used by `ActorQLGYRO`).

Mapped channels: `QeGB→Qe`, `QiGB→Qi`, `GeGB→Γe`, `MtGB→Πi` — the set FUSE flux-matches against
(`GACODE.sources_to_gyrobohm`). The impurity particle flux `GZGB` and turbulent exchange `QieGB`
are not stored.
"""
function write_gacode_fluxes!(actor::ActorFluxMatcherGACODE)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    flux_types = (:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :momentum_flux)

    for (key, model_id, model_name) in ((:turb, :anomalous, "PORTALS (turbulent)"), (:neoc, :neoclassical, "PORTALS (neoclassical)"))
        haskey(actor.flux_data, key) || continue
        jd = actor.flux_data[key]
        fm = jd["fluxes_mean"]
        rho_flux = Float64.(jd["additional_info"]["rho"])
        Qe = Float64.(fm["QeGB"])
        Qi = Float64.(fm["QiGB"])
        Ge = Float64.(fm["GeGB"])
        Mt = Float64.(fm["MtGB"])
        # GB fluxes → GACODE.FluxSolution (PARTICLE_FLUX_i left empty: ion particle flux not mapped)
        sols = [GACODE.FluxSolution{Float64}(Qe[k], Qi[k], Ge[k], Float64[], Mt[k]) for k in eachindex(rho_flux)]

        model = resize!(dd.core_transport.model, model_id; wipe=false)
        model.identifier.name = model_name
        m1d = resize!(model.profiles_1d)
        m1d.grid_flux.rho_tor_norm = rho_flux
        GACODE.flux_gacode_to_imas(flux_types, sols, m1d, eqt, cp1d)
    end

    return actor
end
