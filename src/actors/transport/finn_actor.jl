import TurbulentTransport
import TurbulentTransport: InputTGLF
import GACODE

#= ========= =#
#  ActorFINN  #
#= ========= =#
@actor_parameters_struct ActorFINN{T} begin
    finn_model::Entry{String} = Entry{String}("-", "FINN model filename (BSON)")
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm grid for FINN predictions"; default=0.25:0.1:0.85)
    MXH_modes::Entry{Int} = Entry{Int}("-", "Number of MXH harmonics used when generating InputTGLF geometry for FINN"; default=1)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Warn if inputs are outside FINN training bounds"; default=false)
    evolve_rotation::Switch{Symbol} = Switch{Symbol}(
        [:flux_match, :fixed, :replay],
        "-",
        "Rotation evolution: `:flux_match` (set rotation from predicted VEXB_SHEAR), `:fixed` (leave existing rotation profile alone), `:replay` (core from `replay_dd`, edge from current sim — mirrors `ActorFluxMatcher.evolve_rotation = :replay`)";
        default=:fixed)
    z_max::Entry{Union{T,NamedTuple}} = Entry{Union{T,NamedTuple}}(
        "m⁻¹",
        """
        Maximum allowed normalized gradient (inverse scale length). Can be:
        * Single value: applies to all channels and radii (e.g., 10.0)
        * NamedTuple for spatially varying limits:
          (core=20.0, edge=100.0, rho_transition=0.80)
          Values are constant at 'core' for rho <= rho_transition, then linearly increase to 'edge' at rho=1.0
        """;
        default=10.0,
        check=x -> begin
            if x isa Real
                @assert x > 0.0 "z_max must be positive"
            elseif x isa NamedTuple
                @assert haskey(x, :core) && haskey(x, :edge) && haskey(x, :rho_transition) "Spatially varying z_max must have :core, :edge, :rho_transition"
                @assert x.core > 0.0 && x.edge > 0.0 "core and edge z_max must be positive"
                @assert 0.0 <= x.rho_transition <= 1.0 "rho_transition must be between 0 and 1"
            end
        end
    )
end

mutable struct ActorFINN{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorFINN{P}}
    finn_results::NamedTuple
    actor_replay::ActorReplay{D,P}
end

"""
    ActorFINN(dd::IMAS.dd, act::ParametersAllActors; kw...)

Predicts flux-matched plasma profiles using the FINN (Flux-matcher Inversion Neural Network).

FINN directly predicts the converged TGLF gradients from geometry and source parameters
in a single neural network forward pass (< 1 ms), bypassing the iterative flux-matching loop
that typically takes tens of seconds with TGLF-NN.

Predicted channels (updated from FINN):
- Electron temperature (from RLTS_1)
- Ion temperature (from RLTS_2, applied to all ion species)
- Electron density (from RLNS_1)
- Rotation: see `evolve_rotation` parameter
   - `:flux_match` → set from predicted VEXB_SHEAR
   - `:fixed` → leave existing rotation profile alone
   - `:replay` → core from `replay_dd`, edge from current sim
     (mirrors `ActorFluxMatcher.evolve_rotation == :replay`)

Unpredicted channels (retain experimental/initialized values):
- Individual ion densities (Deuterium, Tritium, impurities)
- Any profiles outside the `rho_transport` grid (edge/pedestal region)
"""
function ActorFINN(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFINN(dd, act.ActorFINN, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFINN(dd::IMAS.dd{D}, par::FUSEparameters__ActorFINN{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorFINN)
    par = OverrideParameters(par; kw...)
    empty_results = (RLTS_1=Float64[], RLTS_2=Float64[], RLNS_1=Float64[], VEXB_SHEAR=Float64[], rho=Float64[])
    # Build a placeholder ActorReplay so the parent struct can be constructed,
    # then reassign with the parent-actor reference so dispatch reaches
    # `_step(::ActorReplay, ::ActorFINN, ::IMAS.dd)` (mirrors ActorFluxMatcher).
    actor_replay = ActorReplay(dd, act.ActorReplay, ActorNoOperation(dd, act.ActorNoOperation))
    actor = ActorFINN(dd, par, empty_results, actor_replay)
    actor.actor_replay = ActorReplay(dd, act.ActorReplay, actor)
    return actor
end

"""
    _step(actor::ActorFINN)

Runs FINN to predict flux-matched gradients and sets plasma profiles.

1. Extracts geometry from equilibrium via InputTGLF
2. Computes sources in gyro-Bohm units from dd.core_sources
3. Runs FINN model to predict TGLF gradients at each rho
4. Converts TGLF gradients to inverse scale lengths (z-profiles)
5. Reconstructs temperature, density profiles from predicted gradients
"""
function _step(actor::ActorFINN{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    # Replay rotation from replay_dd if requested. Runs BEFORE FINN so that
    # the replayed core/edge rotation profile is in place when downstream
    # consumers query `cp1d.rotation_frequency_tor_sonic`. FINN's `_finalize`
    # will then leave rotation untouched (it only writes when `:flux_match`).
    if par.evolve_rotation == :replay
        finalize(step(actor.actor_replay))
    end

    model_filename = endswith(par.finn_model, ".bson") ? par.finn_model[begin:end-5] : par.finn_model
    actor.finn_results = TurbulentTransport.run_finn(
        dd, par.rho_transport;
        model_filename=model_filename,
        warn_nn_train_bounds=par.warn_nn_train_bounds,
        MXH_modes=par.MXH_modes
    )

    return actor
end

"""
    _finalize(actor::ActorFINN)

Converts predicted TGLF gradients to profiles and writes them to dd.core_profiles.

TGLF gradient convention: RLTS_1 = -a/L_Te = -(a/Te) * dTe/dr
Conversion to z-profile: z_Te = (1/Te) * dTe/drho = -RLTS_1 * DRMINDX_LOC

Rotation uses a different conversion because VEXB_SHEAR is not a logarithmic
gradient — it encodes dω/dr with q, rmin, and c_s normalization:
  VEXB_SHEAR = (rmin/q) * (a/c_s) * (-dω/dr)
We solve for dω/drho and integrate inward from the boundary using
`profile_from_rotation_shear_transport`.

Only the channels predicted by FINN are updated:
  - Electron temperature (from RLTS_1)
  - Ion temperature (from RLTS_2)
  - Electron density (from RLNS_1)
  - Rotation (from VEXB_SHEAR, if `evolve_rotation=true`)

Channels NOT predicted by FINN (e.g. individual ion densities like Deuterium)
retain their experimental/initialized values.

Profiles outside the `rho_transport` grid are preserved at their existing values
via `IMAS.profile_from_z_transport`, which only modifies the core region
(ρ ≤ rho_transport[end]) and keeps edge/pedestal profiles intact.
"""
function _finalize(actor::ActorFINN)
    dd = actor.dd
    par = actor.par
    res = actor.finn_results

    cp1d = dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm

    rho_transport = collect(par.rho_transport)

    input_tglfs = InputTGLF(dd, rho_transport, :sat3, true, true; MXH_modes=par.MXH_modes)
    if hasproperty(input_tglfs, :tglfs)
        input_tglfs = input_tglfs.tglfs
    end
    drmindx = [Float64(input_tglfs[i].DRMINDX_LOC) for i in eachindex(rho_transport)]

    # Convert TGLF gradients to z-profiles (inverse scale lengths in rho-space)
    # RLTS_1 = -(a/Te)*dTe/dr  and  z = (1/Te)*dTe/drho = -RLTS * DRMINDX_LOC
    z_Te = -res.RLTS_1 .* drmindx
    z_Ti = -res.RLTS_2 .* drmindx
    z_ne = -res.RLNS_1 .* drmindx

    # Apply z_max clamping
    if par.z_max isa Real
        z_max_val = Float64(par.z_max)
        z_Te = clamp.(z_Te, -z_max_val, z_max_val)
        z_Ti = clamp.(z_Ti, -z_max_val, z_max_val)
        z_ne = clamp.(z_ne, -z_max_val, z_max_val)
    elseif par.z_max isa NamedTuple
        core_limit = par.z_max.core
        edge_limit = par.z_max.edge
        rho_trans = par.z_max.rho_transition
        slope = (edge_limit - core_limit) / (1.0 - rho_trans)
        z_max_profile = [rho_x <= rho_trans ? core_limit : core_limit + slope * (rho_x - rho_trans) for rho_x in rho_transport]
        z_Te = clamp.(z_Te, -z_max_profile, z_max_profile)
        z_Ti = clamp.(z_Ti, -z_max_profile, z_max_profile)
        z_ne = clamp.(z_ne, -z_max_profile, z_max_profile)
    end

    # Electron temperature: update within rho_transport, preserve edge/pedestal
    cp1d.electrons.temperature = IMAS.profile_from_z_transport(
        cp1d.electrons.temperature, rho, rho_transport, z_Te)

    # Ion temperature: update t_i_average within rho_transport, propagate to all ions
    Ti_new = IMAS.profile_from_z_transport(
        cp1d.t_i_average, rho, rho_transport, z_Ti)
    for ion in cp1d.ion
        ion.temperature = Ti_new
    end
    IMAS.unfreeze!(cp1d, :t_i_average)

    # Electron density: update within rho_transport, preserve edge/pedestal
    cp1d.electrons.density_thermal = IMAS.profile_from_z_transport(
        cp1d.electrons.density_thermal, rho, rho_transport, z_ne)
    IMAS.unfreeze!(cp1d.electrons, :density)

    # Ion densities: NOT modified — FINN does not yet predict individual ion density
    # gradients (e.g. Deuterium, Tritium, impurity densities).
    # These retain their experimental/initialized values.

    # Rotation: only update from predicted VEXB_SHEAR when explicitly :flux_match.
    # `:fixed` leaves rotation alone; `:replay` was handled in `_step` via `actor_replay`.
    # From tglf.jl: w0 = -ω, gamma_e = (rmin/q)*dω/dr, VEXB_SHEAR = -gamma_e*(a/c_s)
    # Inverting: dω/drho = VEXB_SHEAR * q / rmin * c_s * DRMINDX
    # where rmin [cm] = a_cm * RMIN_LOC, and c_s [cm/s] from GACODE
    if par.evolve_rotation == :flux_match && !isempty(res.VEXB_SHEAR)
        eqt1d = dd.equilibrium.time_slice[].profiles_1d
        rmin_full = GACODE.r_min_core_profiles(eqt1d, rho)
        a_cm = rmin_full[end]

        cp_gridpoints = [IMAS.argmin_abs(rho, rho_x) for rho_x in rho_transport]
        Te_transport = cp1d.electrons.temperature[cp_gridpoints]
        c_s = GACODE.c_s.(Te_transport)
        q_loc = [Float64(input_tglfs[i].Q_LOC) for i in eachindex(rho_transport)]
        rmin_loc = [Float64(input_tglfs[i].RMIN_LOC) for i in eachindex(rho_transport)]

        dw_drho = @. res.VEXB_SHEAR * q_loc / (a_cm * rmin_loc) * c_s * drmindx

        cp1d.rotation_frequency_tor_sonic = IMAS.profile_from_rotation_shear_transport(
            cp1d.rotation_frequency_tor_sonic, rho, rho_transport, dw_drho)
    end

    IMAS.intrinsic_sources!(dd)

    return actor
end

"""
    _step(replay_actor::ActorReplay, actor::ActorFINN, replay_dd::IMAS.dd)

Replay rotation profile from `replay_dd` into the current `dd` when
`actor.par.evolve_rotation == :replay`. Mirrors the rotation block of
`ActorFluxMatcher`'s replay step (core from replay, edge from current sim,
shifted at `rho_transport[end]` for continuity). Both `rotation_frequency_tor_sonic`
and per-ion `rotation_frequency_tor` are copied — the latter is the measured
quantity, and copying only sonic rotation would let it get recomputed from
simulated pressure gradients via the IMAS expression.
"""
function _step(replay_actor::ActorReplay, actor::ActorFINN, replay_dd::IMAS.dd)
    par = actor.par
    par.evolve_rotation == :replay || return replay_actor

    dd = actor.dd
    time0 = dd.global_time
    cp1d = dd.core_profiles.profiles_1d[time0]
    replay_cp1d = replay_dd.core_profiles.profiles_1d[time0]
    rho = cp1d.grid.rho_tor_norm

    # Use the outer edge of the FINN transport grid as the blend point — this is
    # the analog of FluxMatcher's `rho_nml`.
    rho_blend = collect(par.rho_transport)[end]
    i_blend = IMAS.argmin_abs(rho, rho_blend)

    # Sonic rotation: core from replay, edge from simulation, shift to match at i_blend
    ω_core = deepcopy(replay_cp1d.rotation_frequency_tor_sonic)
    ω_edge = IMAS.freeze!(cp1d, :rotation_frequency_tor_sonic)
    ω_core[i_blend+1:end] = ω_edge[i_blend+1:end]
    ω_core[1:i_blend] = ω_core[1:i_blend] .- ω_core[i_blend] .+ ω_edge[i_blend]
    cp1d.rotation_frequency_tor_sonic = ω_core

    # Ion rotation: same blending (core from replay, edge from simulation)
    for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
        if IMAS.hasdata(replay_ion, :rotation_frequency_tor)
            ω_ion_core = deepcopy(replay_ion.rotation_frequency_tor)
            ω_ion_edge = IMAS.freeze!(ion, :rotation_frequency_tor)
            ω_ion_core[i_blend+1:end] = ω_ion_edge[i_blend+1:end]
            ω_ion_core[1:i_blend] = ω_ion_core[1:i_blend] .- ω_ion_core[i_blend] .+ ω_ion_edge[i_blend]
            ion.rotation_frequency_tor = ω_ion_core
        end
    end

    return replay_actor
end
