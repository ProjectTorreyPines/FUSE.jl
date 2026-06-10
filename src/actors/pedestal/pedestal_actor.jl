#= ============= =#
#  ActorPedestal  #
#= ============= =#
@actor_parameters_struct ActorPedestal{T} begin
    #== common pedestal parameters==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
    T_ratio_pedestal::Entry{T} =
        Entry{T}("-", "Ratio of ion to electron temperatures (or rho at which to sample for that ratio, if negative; or rho_nml-(rho_ped-rho_nml) if 0.0)"; default=1.0)
    Te_sep::Entry{T} = Entry{T}("-", "Separatrix electron temperature"; default=80.0, check=x -> @assert x > 0 "Te_sep must be > 0")
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    ne_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_from::Switch{Symbol} = switch_get_from(:zeff_ped)
    mode_transitions::Entry{Dict{Float64,Symbol}} = Entry{Dict{Float64,Symbol}}(
        "s",
        "Times at which the plasma transitions to a given mode [:L_mode, :H_mode]. If missing, the L-H transition will be based on `IMAS.satisfies_h_mode_conditions(dd)`."
    )
    #== actor parameters==#
    density_match::Switch{Symbol} = Switch{Symbol}([:ne_line, :ne_ped], "-", "Matching density based on ne_ped or line averaged density"; default=:ne_ped)
    model::Switch{Symbol} = Switch{Symbol}([:EPED, :WPED, :dynamic, :analytic, :replay, :none], "-", "Pressure edge model"; default=:EPED)
    rotation_model::Switch{Symbol} = Switch{Symbol}(
        [:linear, :replay, :nn_pedestal, :none],
        "-",
        "Rotation edge model: `:linear` (linear edge), `:replay` (core+edge from replay_dd), `:nn_pedestal` (use NN-predicted ω_ped, requires `ne_from=:nn_predictor`), `:none` (do nothing — preserve existing rotation profile)";
        default=:none)
    #== L to H and H to L transition model ==#
    tau_t::Entry{T} = Entry{T}("s", "Edge temperature LH transition tanh evolution time (95% of full transition)")
    tau_n::Entry{T} = Entry{T}("s", "Edge density LH transition tanh evolution time (95% of full transition)")
    density_ratio_L_over_H::Entry{T} = Entry{T}("-", "n_Lmode / n_Hmode")
    zeff_ratio_L_over_H::Entry{T} = Entry{T}("-", "zeff_Lmode / zeff_Hmode")
    #== nn_predictor FPE source ==#
    fpe_source::Switch{Symbol} = Switch{Symbol}([:zmq, :dd], "-", "Source of FPE actuator inputs for nn_predictor: `:zmq` (from ActorZMQ via dd._aux) or `:dd` (directly from dd fields, no ZMQ required)"; default=:zmq)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorPedestal{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorPedestal{P}}
    act::ParametersAllActors{P}
    ped_actor::Union{ActorWPED{D,P},ActorEPED{D,P},ActorAnalyticPedestal{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
    wped_actor::ActorWPED{D,P}
    eped_actor::ActorEPED{D,P}
    analytic_actor::ActorAnalyticPedestal{D,P}
    replay_actor::Union{ActorReplay{D,P},ActorNoOperation{D,P}}
    noop_actor::ActorNoOperation{D,P}
    state::Vector{Symbol}
    t_lh::Float64
    t_hl::Float64
    previous_time::Float64
    cp1d_transition::IMAS.core_profiles__profiles_1d{D}
    nn_predictor::Union{Nothing,PedestalNN}
    nn_prediction::Union{Nothing,NamedTuple}
end

"""
    ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)

Comprehensive pedestal modeling with support for multiple models and L-H mode transitions.

This compound actor manages pedestal physics by selecting from available pedestal models
and handling mode transitions between L-mode and H-mode operation. It coordinates multiple
specialized pedestal actors and provides dynamic transition capabilities.

Available pedestal models:
- `:EPED`: EPED neural network model for pedestal predictions
- `:WPED`: Width-based energy balance pedestal model
- `:analytic`: Analytic scaling laws for spherical tokamaks
- `:dynamic`: Time-dependent L-H transitions with smoothing
- `:replay`: Replays pedestal data from experimental reference
- `:none`: No pedestal modifications

Key features:
- Automatic L-H mode detection based on power threshold criteria
- Dynamic mode transitions with configurable time constants
- Density matching options (pedestal or line-averaged)
- Rotation profile modeling (linear, experimental replay)
- Consistent Ti/Te ratio handling across all models

Mode transition physics:
- Supports user-defined transition times or automatic power threshold detection
- Smooth temporal evolution using tanh functions with configurable time scales
- Separate evolution times for temperature and density transitions
- Configurable density and Zeff ratios between L-mode and H-mode
"""
function ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPedestal(dd, act.ActorPedestal, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorPedestal(dd::IMAS.dd{D}, par::FUSEparameters__ActorPedestal{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorPedestal)
    par = OverrideParameters(par; kw...)
    eped_actor =
        ActorEPED(dd, act.ActorEPED; par.rho_nml, par.rho_ped, par.T_ratio_pedestal, par.Te_sep, par.ip_from, par.βn_from, ne_from=:core_profiles, zeff_from=:core_profiles)
    wped_actor =
        ActorWPED(dd, act.ActorWPED; par.rho_nml, par.rho_ped, par.T_ratio_pedestal, par.Te_sep, par.ip_from, par.βn_from, ne_from=:core_profiles, zeff_from=:core_profiles)
    analytic_actor =
        ActorAnalyticPedestal(
            dd,
            act.ActorAnalyticPedestal;
            par.rho_nml,
            par.rho_ped,
            par.T_ratio_pedestal,
            par.Te_sep,
            par.ip_from,
            par.βn_from,
            ne_from=:core_profiles,
            zeff_from=:core_profiles
        )
    noop = ActorNoOperation(dd, act.ActorNoOperation)
    actor = ActorPedestal(dd, par, act, noop, wped_actor, eped_actor, analytic_actor, noop, noop, Symbol[], -Inf, -Inf, -Inf, IMAS.core_profiles__profiles_1d{D}(), nothing, nothing)
    actor.replay_actor = ActorReplay(dd, act.ActorReplay, actor)
    return actor
end

"""
    _step(actor::ActorPedestal)

Orchestrates pedestal model selection and mode transition logic.

The step function manages the complex workflow of:
1. Determining current plasma mode (L-mode, H-mode) from power balance or user input
2. Selecting and running the appropriate pedestal model
3. Handling dynamic L-H transitions with proper temporal smoothing
4. Applying rotation models if requested
5. Updating plasma profiles with pedestal boundary conditions

For dynamic mode transitions, tracks transition times and applies gradual profile
evolution to avoid numerical discontinuities.
"""
function _step(actor::ActorPedestal{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    # Run NN predictor once per time step (provides ne_ped, Te_ped, Ti_ped,
    # T_rot_ped, and the L/H classifier in a single ensemble call). Live ZMQ
    # actuator signals (stored in dd._aux by ActorZMQ.receive!) are mapped onto
    # the FPE channels; channels with no live source stay at the training mean
    # via the per-channel z-score in normalized_signals.
    if par.ne_from == :nn_predictor
        if actor.nn_predictor === nothing
            actor.nn_predictor = load_pedestal_nn()
            @info "ActorPedestal: loaded NN pedestal predictor ensemble ($(length(actor.nn_predictor.bundles)) bundles) from $(actor.nn_predictor.onnx_dir)"
        end
        sequences, signal_mask = build_fpe_sequences_from_aux(actor.nn_predictor, dd; use_dd=(par.fpe_source == :dd))
        # Live MSE history (from `dd._aux[:nn_history_buffer]` populated by
        # end-of-shot `push_shot_history!` calls). When the buffer is empty we
        # silently fall back to mean_normalized_history, equivalent to the
        # pre-buffer behavior. norm_params_path stays `nothing` for now until
        # the per-block source registries (compute_shot_stats) are wired —
        # raw stats are biased but the buffer plumbing is exercised end-to-end.
        history = mse_history_from_aux(dd, actor.nn_predictor)
        actor.nn_prediction = predict_pedestal(actor.nn_predictor; sequences, signal_mask, history)
    end

    if !ismissing(par, :mode_transitions)
        causal_transition_time = IMAS.nearest_causal_time(sort!(collect(keys(par.mode_transitions))), dd.global_time).causal_time
        mode = par.mode_transitions[causal_transition_time]
    elseif par.ne_from == :nn_predictor && actor.nn_prediction !== nothing
        mode = actor.nn_prediction.is_h_mode ? :H_mode : :L_mode
    elseif IMAS.satisfies_h_mode_conditions(dd; threshold_multiplier=1.2)
        mode = :H_mode
    elseif !IMAS.satisfies_h_mode_conditions(dd; threshold_multiplier=0.8)
        mode = :L_mode
    elseif isempty(actor.state)
        if IMAS.satisfies_h_mode_conditions(dd)
            mode = :H_mode
        else
            mode = :L_mode
        end
    else
        mode = actor.state[end]
    end
    push!(actor.state, mode)
    @ddtime(dd.summary.global_quantities.h_mode.value = Int(mode == :H_mode))

    if par.model == :none
        actor.ped_actor = actor.noop_actor
        finalize(step(actor.ped_actor))

    elseif par.model == :replay
        actor.ped_actor = actor.replay_actor
        finalize(step(actor.ped_actor))

    else
        if par.model == :EPED
            actor.ped_actor = actor.eped_actor
            run_selected_pedestal_model(actor; density_factor=1.0, zeff_factor=1.0)

        elseif par.model == :WPED
            actor.ped_actor = actor.wped_actor
            run_selected_pedestal_model(actor; density_factor=1.0, zeff_factor=1.0)

        elseif par.model == :analytic
            actor.ped_actor = actor.analytic_actor
            run_selected_pedestal_model(actor; density_factor=1.0, zeff_factor=1.0)

        elseif par.model == :dynamic
            @assert par.ne_from in (:pulse_schedule, :nn_predictor) ":dynamic pedestal model requires `act.ActorPedestal.ne_from` ∈ (:pulse_schedule, :nn_predictor)"
            @assert actor.previous_time < dd.global_time "subsequent calls to :dynamic pedestal model require dd.global_time advance"

            if length(actor.state) < 2
                # initialization
                actor.t_lh = -Inf
                actor.t_hl = -Inf
                actor.cp1d_transition = deepcopy(cp1d)

            elseif length(actor.state) >= 2 && actor.state[end-1:end] == [:L_mode, :H_mode]
                # L to H
                actor.t_lh = (actor.previous_time + dd.global_time) / 2.0
                actor.cp1d_transition = deepcopy(cp1d)

            elseif length(actor.state) >= 2 && actor.state[end-1:end] == [:H_mode, :L_mode]
                # H to L
                actor.t_hl = (actor.previous_time + dd.global_time) / 2.0
                actor.cp1d_transition = deepcopy(cp1d)
            end

            if mode == :L_mode
                # L-mode
                α_t = LH_dynamics(par.tau_t, actor.t_hl, dd.global_time) # from 0 -> 1
                α_n = LH_dynamics(par.tau_n, actor.t_hl, dd.global_time) # from 0 -> 1

                actor.ped_actor = actor.wped_actor
                density_factor = 1.0 * (1 - α_n) + par.density_ratio_L_over_H * α_n
                zeff_factor = 1.0 * (1 - α_n) + par.zeff_ratio_L_over_H * α_n

                run_selected_pedestal_model(actor; density_factor, zeff_factor)

                Te_now = (1 .- α_t) .* actor.cp1d_transition.electrons.temperature .+ α_t .* cp1d.electrons.temperature
                Ti_now = (1 .- α_t) .* actor.cp1d_transition.ion[1].temperature .+ α_t .* cp1d.ion[1].temperature

                cp1d.electrons.temperature = Te_now
                for ion in cp1d.ion
                    ion.temperature = Ti_now
                end

            else
                # H-mode
                α_t = LH_dynamics(par.tau_t, actor.t_lh, dd.global_time) # from 0 -> 1
                α_n = LH_dynamics(par.tau_n, actor.t_lh, dd.global_time) # from 0 -> 1

                actor.ped_actor = actor.eped_actor
                density_factor = par.density_ratio_L_over_H * (1 - α_n) + 1.0 * α_n
                zeff_factor = par.zeff_ratio_L_over_H * (1 - α_n) + 1.0 * α_n

                run_selected_pedestal_model(actor; density_factor, zeff_factor)

                Te_now = (1 .- α_t) .* actor.cp1d_transition.electrons.temperature .+ α_t .* cp1d.electrons.temperature
                Ti_now = (1 .- α_t) .* actor.cp1d_transition.ion[1].temperature .+ α_t .* cp1d.ion[1].temperature

                cp1d.electrons.temperature = Te_now
                for ion in cp1d.ion
                    ion.temperature = Ti_now
                end
            end

        end

        if par.rotation_model == :linear
            # linear pedestal rotation with zero boundary condition at the edge
            rho = cp1d.grid.rho_tor_norm
            i_nml = IMAS.argmin_abs(rho, par.rho_nml)
            i_ped = IMAS.argmin_abs(rho, par.rho_ped)
            ω_core = IMAS.freeze!(cp1d.ion[1], :rotation_frequency_tor)
            if i_nml == i_ped
                dωdr_nml = IMAS.gradient(rho, -ω_core; method=:backward)[i_nml]
            else
                dωdr_nml = (ω_core[i_nml] - ω_core[i_ped]) / (rho[i_ped] - rho[i_nml])
            end
            ω_edge_linear = (1.0 .- rho) * dωdr_nml
            ω_core[i_nml+1:end] = ω_edge_linear[i_nml+1:end]
            ω_core[1:i_nml] = ω_core[1:i_nml] .- ω_core[i_nml] .+ ω_edge_linear[i_nml]
            for ion in cp1d.ion
                ion.rotation_frequency_tor = ω_core
            end
            IMAS.ωtor2sonic!(cp1d)

        elseif par.rotation_model == :replay
            # Edge from replay, core from simulation (opposite of FluxMatcher)
            # NOTE: We must also copy ion.rotation_frequency_tor, not just sonic rotation,
            # because ion rotation is the measured quantity.
            time0 = dd.global_time
            rho = cp1d.grid.rho_tor_norm
            replay_cp1d = actor.replay_actor.replay_dd.core_profiles.profiles_1d[time0]
            i_nml = IMAS.argmin_abs(rho, par.rho_nml)
            i_ped = IMAS.argmin_abs(rho, par.rho_ped)

            # Sonic rotation: core from simulation, edge from replay, shift to match
            ω_core = IMAS.freeze!(cp1d, :rotation_frequency_tor_sonic)
            ω_edge = replay_cp1d.rotation_frequency_tor_sonic
            ω_core[i_nml+1:end] = ω_edge[i_nml+1:end]
            ω_core[1:i_nml] = ω_core[1:i_nml] .- ω_core[i_nml] .+ ω_edge[i_nml]
            cp1d.rotation_frequency_tor_sonic = ω_core

            # Ion rotation: same blending (core from simulation, edge from replay)
            for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
                if IMAS.hasdata(replay_ion, :rotation_frequency_tor)
                    ω_ion_core = IMAS.freeze!(ion, :rotation_frequency_tor)
                    ω_ion_edge = replay_ion.rotation_frequency_tor
                    ω_ion_core[i_nml+1:end] = ω_ion_edge[i_nml+1:end]
                    ω_ion_core[1:i_nml] = ω_ion_core[1:i_nml] .- ω_ion_core[i_nml] .+ ω_ion_edge[i_nml]
                    ion.rotation_frequency_tor = ω_ion_core
                end
            end
        end

    end

    actor.previous_time = dd.global_time

    return actor
end

function _finalize(actor::ActorPedestal{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd

    cp1d = dd.core_profiles.profiles_1d[]
    summary_ped = dd.summary.local.pedestal
    rho = cp1d.grid.rho_tor_norm

    IMAS.enforce_quasi_neutrality!(cp1d, :electrons)

    position = 1 - IMAS.pedestal_tanh_width_half_maximum(rho, cp1d.electrons.temperature)
    @ddtime summary_ped.position.rho_tor_norm = position
    @ddtime summary_ped.n_e.value = IMAS.interp1d(rho, cp1d.electrons.density_thermal).(position)
    @ddtime summary_ped.zeff.value = IMAS.interp1d(rho, cp1d.zeff).(position)
    @ddtime summary_ped.t_e.value = IMAS.interp1d(rho, cp1d.electrons.temperature).(position)
    @ddtime summary_ped.t_i_average.value = IMAS.interp1d(rho, cp1d.t_i_average).(position)

    return actor
end

"""
    LH_dynamics(τ::Float64, t_LH::Float64, t_now::Float64)

Returns a parameter that follows a tanh like response where τ represent the value of 0.95 @ τ time starting from t_LH
"""
function LH_dynamics(τ::Float64, t_LH::Float64, t_now::Float64)
    if t_LH == -Inf
        return 1.0
    elseif t_now <= t_LH
        return 0.0
    end
    α = tanh.((2pi .* (t_now .- t_LH .- τ / 2.0)) ./ τ) / 2.0 + 0.5
    α0 = tanh.((2pi .* (.-τ / 2.0)) ./ τ) / 2.0 + 0.5
    α = (α .- α0) ./ (1 - α0)
    return α
end

"""
    pedestal_density_tanh(dd::IMAS.dd, par::OverrideParameters{P,FUSEparameters__ActorPedestal{P}}; density_factor::Float64, zeff_factor::Float64) where {P<:Real}

The edge density must be defined independently of the pedestal model

The EPED and WPED models only operate on the temperature profiles
"""
function pedestal_density_tanh(dd::IMAS.dd, par::OverrideParameters{P,FUSEparameters__ActorPedestal{P}};
                               density_factor::Float64, zeff_factor::Float64,
                               nn_prediction::Union{Nothing,NamedTuple}=nothing) where {P<:Real}
    cp1d = dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm

    # Throughout FUSE, the "pedestal" values are defined at rho=0.9
    rho09 = 0.9

    # density pedestal width to match the existing temperature pedestal width
    w_ped = IMAS.pedestal_tanh_width_half_maximum(rho, cp1d.electrons.temperature)

    ne_old = copy(cp1d.electrons.density_thermal)
    if par.ne_from == :nn_predictor
        @assert nn_prediction !== nothing "ne_from=:nn_predictor requires nn_prediction (call predict_density first)"
        ne_ped_1e19 = sum(nn_prediction.predictions_physical) / length(nn_prediction.predictions_physical)
        ne_ped = Float64(ne_ped_1e19) * 1e19 * density_factor
        @info "ActorPedestal: nn_predictor ne_ped = $(round(ne_ped_1e19; digits=3)) x 10^19 m^-3"
    else
        ne_ped = IMAS.get_from(dd, Val(:ne_ped), par.ne_from, rho09) * density_factor
    end
    cp1d.electrons.density_thermal[end] = ne_ped / 4.0
    ne = IMAS.blend_core_edge_Hmode(cp1d.electrons.density_thermal, rho, ne_ped, w_ped, par.rho_nml, par.rho_ped; method=:scale)
    cp1d.electrons.density_thermal = ne = IMAS.ped_height_at_09(rho, ne, ne_ped)
    IMAS.unfreeze!(cp1d.electrons, :density)
    ratio = ne ./ ne_old

    for ion in cp1d.ion
        if !ismissing(ion, :density_thermal)
            ion.density_thermal = ion.density_thermal .* ratio
            ni_ped = IMAS.interp1d(rho, ion.density_thermal).(rho09)
            ion.density_thermal[end] = ni_ped / 4.0
            ni = IMAS.blend_core_edge_Hmode(ion.density_thermal, rho, ni_ped, w_ped, par.rho_nml, par.rho_ped; method=:scale)
            ion.density_thermal = IMAS.ped_height_at_09(rho, ni, ni_ped)
            IMAS.unfreeze!(ion, :density)
        end
    end

    #NOTE: Zeff can change after a pedestal actor is run, even though actors like EPED and WPED only operate on the temperature profiles.
    # This is because in FUSE the calculation of Zeff is temperature dependent.
    zeff_ped = IMAS.get_from(dd, Val(:zeff_ped), par.zeff_from, rho09) * zeff_factor
    IMAS.scale_ion_densities_to_target_zeff!(cp1d, rho09, zeff_ped)

    return nothing
end

"""
    build_fpe_sequences_from_aux(nn::PedestalNN, dd::IMAS.dd; T::Integer=200) -> (Matrix{Float32}, Vector{Float32})

Build a `(T, 32)` raw-physical FPE input matrix and a length-32 `signal_mask`
from the live ZMQ signals stored in `dd._aux` (populated by `ActorZMQ.receive!`)
and standard IMAS fields. Channels that have no live data source are filled
with the per-channel training mean (so they z-score to zero inside
[`predict_pedestal`](@ref), equivalent to the "missing" / mean-input fallback)
and their `signal_mask` bit stays 0.

DIII-D ZMQ → FPE channel mapping:
- `dd._aux[:zmq_Ip_avg].values`           -> `ip`         (plasma current, A)
- `dd._aux[:zmq_pr15v].values`            -> `ipspr15v`   (P-coil 15 V regulator current, A)
- `dd._aux[:zmq_Pohm].values`             -> `pohm`       (ohmic power, W; passthrough)
- `dd._aux[:zmq_Pnbi].values`             -> `pinj`       (total NBI power, W → MW for FPE)
- `dd._aux[:zmq_gas[a-e]_cal].values`     -> `gas[a-e]_cal` (calibrated gas flow)
- `dd._aux[:zmq_I_coil].values[1]`        -> `ecoila`     (E-coil A bank, PCECOILA)
- `dd._aux[:zmq_I_coil].values[4]`        -> `ecoilb`     (E-coil B bank, PCECOILB)
- `dd._aux[:zmq_I_coil].values[7..15]`    -> `f1a..f9a`   (F-coil A bank currents, A)
- `dd._aux[:zmq_I_coil].values[16..24]`   -> `f1b..f9b`   (F-coil B bank currents, A)
- `dd.equilibrium.vacuum_toroidal_field.b0[end]` -> `bt` (T)

`I_coil` is a 24-element vector of DIII-D PCS coil-current pointnames in this
order (confirmed against the `PCSpcsRtnetCoilNames` MATLAB lookup):

| idx | pointname | FPE channel |
|----:|-----------|-------------|
|  1  | PCECOILA  | `ecoila`    |
|  2  | PCE89DN   | (unused)    |
|  3  | PCE567UP  | (unused)    |
|  4  | PCECOILB  | `ecoilb`    |
|  5  | PCE89UP   | (unused)    |
|  6  | PCE567DN  | (unused)    |
|  7  | PCF1A     | `f1a`       |
| ... | ...       | ...         |
| 15  | PCF9A     | `f9a`       |
| 16  | PCF1B     | `f1b`       |
| ... | ...       | ...         |
| 24  | PCF9B     | `f9b`       |

PedestalPredictor's FPE only uses the two E-coils (`ecoila`, `ecoilb`); the
internal C-coil segments at indices 2,3,5,6 (`PCE89DN`, `PCE567UP`, `PCE89UP`,
`PCE567DN`) are not consumed. Bounds checks (`length(v) >= idx`) make the
mapping safe against shorter `I_coil` vectors.

All FPE channels are now wired to live sources.

Returns:
- `sequences::Matrix{Float32}` — `(T, 32)`, raw physical units, broadcast across time.
- `signal_mask::Vector{Float32}` — length 32, 1.0 where a live value was found, 0.0 otherwise.
"""
function build_fpe_sequences_from_aux(nn::PedestalNN, dd::IMAS.dd; T::Integer=200, use_dd::Bool=false)
    # Initialize every channel to its training mean so an un-mapped channel
    # z-scores to exactly zero (matching mean_normalized_history's convention).
    sequences = repeat(reshape(nn.signal_means, 1, 32), T, 1)
    mask = zeros(Float32, 32)

    aux = getfield(dd, :_aux)
    t_now = dd.global_time
    mapped = String[]

    function _set_channel!(name::AbstractString, value::Real)
        idx = fpe_signal_index(nn, name)
        if idx == 0
            @warn "build_fpe_sequences_from_aux: FPE channel \"$name\" not in nn.signal_names; skipping"
            return
        end
        sequences[:, idx] .= Float32(value)
        mask[idx] = 1f0
        push!(mapped, name)
        return
    end

    # Latest causal sample (t_i <= t_now); fall back to the first sample if all
    # entries are in the future (e.g. just after a time-rewind).
    function _aux_value_at(key::Symbol)
        haskey(aux, key) || return nothing
        rec = aux[key]
        (hasproperty(rec, :times) && hasproperty(rec, :values)) || return nothing
        isempty(rec.times) && return nothing
        idx = findlast(τ -> τ <= t_now, rec.times)
        idx === nothing && (idx = 1)
        return rec.values[idx]
    end

    # Scalar mappings.
    # Units reconciliation between the ZMQ wire (A / W / N·m) and the FPE
    # training set (inspect per-bundle `normalization_params.json::means/stds`
    # alongside `onnx_models/fpe_pre_normalization_params.json`):
    # - `ip`       : wire A, training A (edensfit89 bundle mean≈3.5e4, std≈2e5)  -> passthrough
    # - `ipspr15v` : wire A, training MA (bundle mean≈0.07, std≈1.06)            -> /1e6
    # - `pohm`     : wire W, training W (bundle mean≈4e5, std≈4e5)               -> passthrough
    # - `pinj`     : wire W, training MW (bundle mean≈2.8, std≈194)              -> /1e6
    # - `ech_total`: wire W, training MW (bundle mean≈2.2, std≈6.0)              -> /1e6
    # - `tinj`     : wire N·m, training N·m (bundle mean≈2.2, std≈2.4)           -> passthrough
    if use_dd
        # dd path: read FPE inputs directly from dd fields (no ZMQ required)
        t_now = dd.global_time

        # ip
        if !isempty(dd.equilibrium.time_slice)
            _set_channel!("ip", dd.equilibrium.time_slice[].global_quantities.ip)
        end

        # pohm — from core_sources ohmic source
        ohmic_srcs = IMAS.findall(dd.core_sources.source, "identifier.name" => "ohmic")
        if !isempty(ohmic_srcs) && !isempty(ohmic_srcs[1].profiles_1d)
            _set_channel!("pohm", IMAS.total_power_source(ohmic_srcs[1].profiles_1d[end]))
        end

        # pinj — total NBI power from pulse_schedule
        if !isempty(dd.pulse_schedule.nbi.unit)
            Pnbi = 0.0
            for unit in dd.pulse_schedule.nbi.unit
                if !ismissing(unit.power, :reference) && !isempty(unit.power.reference.time)
                    Pnbi += IMAS.interp1d(unit.power.reference.time, unit.power.reference.data, :constant)(t_now)
                end
            end
            _set_channel!("pinj", Pnbi / 1e6)
        end

        # pech — total EC power from ec_launchers
        if !isempty(dd.ec_launchers.beam)
            Pech = sum(
                IMAS.interp1d(beam.power_launched.time, beam.power_launched.data, :constant)(t_now)
                for beam in dd.ec_launchers.beam
            )
            _set_channel!("ech_total", Pech / 1e6)
        end

        # ipspr15v, gasa..gase_cal, ecoila, ecoilb, f1a..f9b — not available from dd; fall back to training mean

    else
        # zmq path: read FPE inputs from dd._aux populated by ActorZMQ.receive!
        let v = _aux_value_at(:zmq_Ip_avg);  v === nothing || _set_channel!("ip",       v); end
        let v = _aux_value_at(:zmq_pr15v);   v === nothing || _set_channel!("ipspr15v", v / 1e6); end
        let v = _aux_value_at(:zmq_Pohm);    v === nothing || _set_channel!("pohm", v); end
        let v = _aux_value_at(:zmq_Pnbi);    v === nothing || _set_channel!("pinj", v / 1e6); end
        let v = _aux_value_at(:zmq_Pech);    v === nothing || _set_channel!("ech_total", v / 1e6); end
        for (k, name) in zip(
                (:zmq_gasa_cal, :zmq_gasb_cal, :zmq_gasc_cal, :zmq_gasd_cal, :zmq_gase_cal),
                ("gasa_cal",    "gasb_cal",    "gasc_cal",    "gasd_cal",    "gase_cal"))
            v = _aux_value_at(k)
            v === nothing || _set_channel!(name, v)
        end

        # PCS coil currents
        let v = _aux_value_at(:zmq_I_coil)
            if v !== nothing
                length(v) >= 1 && _set_channel!("ecoila", v[1])
                length(v) >= 4 && _set_channel!("ecoilb", v[4])
                f_names = ("f1a","f2a","f3a","f4a","f5a","f6a","f7a","f8a","f9a",
                           "f1b","f2b","f3b","f4b","f5b","f6b","f7b","f8b","f9b")
                for (k, name) in enumerate(f_names)
                    length(v) >= 6 + k || break
                    _set_channel!(name, v[6+k])
                end
            end
        end
    end

    # tinj — NBI torque from core_sources for both paths (NBI identifier index = 2)
    tinj = 0.0
    for src in dd.core_sources.source
        if src.identifier.index == 2 && !isempty(src.global_quantities)
            gq = src.global_quantities[]
            if !ismissing(gq, :torque_tor)
                tinj += gq.torque_tor
            end
        end
    end
    _set_channel!("tinj", tinj)

    # bt — from dd.equilibrium for both paths
    if !isempty(dd.equilibrium.vacuum_toroidal_field.b0)
        _set_channel!("bt", dd.equilibrium.vacuum_toroidal_field.b0[end])
    end

    if isempty(mapped)
        @debug "ActorPedestal: nn_predictor — no live FPE channels available; running on training means"
    else
        @debug "ActorPedestal: nn_predictor mapped $(length(mapped)) live FPE channels: $(mapped)"
    end

    return sequences, mask
end

"""
    pedestal_nn_apply!(actor::ActorPedestal)

When `par.ne_from == :nn_predictor` and an NN ensemble prediction is attached
to the actor, blend the NN-predicted pedestal temperatures and rotation into
`cp1d` so that downstream consumers (and `_finalize`'s `dd.summary.local.pedestal`
writes) reflect them.

- `nn_prediction.te_ped`     (keV)   -> `cp1d.electrons.temperature`
- `nn_prediction.ti_ped`     (keV)   -> every `cp1d.ion[*].temperature`
- `nn_prediction.t_rot_ped`  (krad/s) -> every `cp1d.ion[*].rotation_frequency_tor`,
   ONLY when `par.rotation_model == :nn_pedestal` (explicit opt-in). The default
   `:none` is a true no-op — rotation is left untouched, mirroring `ActorFluxMatcher`'s
   `evolve_rotation == :fixed` semantics. Use `:linear`/`:replay` for the historical
   rotation BCs, or `:nn_pedestal` to drive the rotation pedestal off the NN.

Each NN field is reduced to a scalar pedestal value by averaging over the
FPE window — the same convention already used for `predictions_physical`
in [`pedestal_density_tanh`](@ref). NaN-filled fields (e.g. when a bundle
was not loaded) are skipped silently.
"""
function pedestal_nn_apply!(actor::ActorPedestal)
    par = actor.par
    nn = actor.nn_prediction
    (par.ne_from == :nn_predictor && nn !== nothing) || return actor

    cp1d = actor.dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm
    rho09 = 0.9
    w_ped = IMAS.pedestal_tanh_width_half_maximum(rho, cp1d.electrons.temperature)

    # Snapshot Ti/Te ratio *before* we mutate Te so the Ti separatrix boundary
    # stays consistent with the post-ped_actor profiles.
    Ti_over_Te = ti_te_ratio(cp1d, par.T_ratio_pedestal, par.rho_nml, par.rho_ped)

    function _trace_mean(x)
        v = filter(isfinite, x)
        return isempty(v) ? NaN : sum(v) / length(v)
    end

    Te_ped_keV = _trace_mean(nn.te_ped)
    if isfinite(Te_ped_keV)
        Te_ped_eV = Te_ped_keV * 1e3
        Te = copy(cp1d.electrons.temperature)
        Te[end] = par.Te_sep
        Te = IMAS.blend_core_edge_Hmode(Te, rho, Te_ped_eV, w_ped, par.rho_nml, par.rho_ped; method=:scale)
        cp1d.electrons.temperature = IMAS.ped_height_at_09(rho, Te, Te_ped_eV)
        @info "ActorPedestal: nn_predictor Te_ped = $(round(Te_ped_keV; digits=3)) keV"
    end

    Ti_ped_keV = _trace_mean(nn.ti_ped)
    if isfinite(Ti_ped_keV)
        Ti_ped_eV = Ti_ped_keV * 1e3
        Ti_sep = par.Te_sep * Ti_over_Te
        for ion in cp1d.ion
            if !ismissing(ion, :temperature)
                Ti = copy(ion.temperature)
                Ti[end] = Ti_sep
                Ti = IMAS.blend_core_edge_Hmode(Ti, rho, Ti_ped_eV, w_ped, par.rho_nml, par.rho_ped; method=:scale)
                ion.temperature = IMAS.ped_height_at_09(rho, Ti, Ti_ped_eV)
            end
        end
        @info "ActorPedestal: nn_predictor Ti_ped = $(round(Ti_ped_keV; digits=3)) keV"
    end

    # rotation_frequency_tor can be negative, so we cannot use blend_core_edge_Hmode
    # (which uses log internally — see the :replay rotation branch in _step).
    # Use a simple scale-to-target at rho=0.9, falling back to a constant offset
    # if the existing trace passes through zero at the pedestal foot.
    # We also refresh cp1d.rotation_frequency_tor_sonic so FINN's
    # `profile_from_rotation_shear_transport` sees the NN pedestal value as BC.
    T_rot_krads = _trace_mean(nn.t_rot_ped)
    if isfinite(T_rot_krads) && par.rotation_model == :nn_pedestal
        ω_ped = T_rot_krads * 1e3
        any_ion_rot = false
        for ion in cp1d.ion
            if !ismissing(ion, :rotation_frequency_tor)
                ω = copy(ion.rotation_frequency_tor)
                ω_at_09 = IMAS.interp1d(rho, ω).(rho09)
                ω = iszero(ω_at_09) ? ω .+ ω_ped : ω .* (ω_ped / ω_at_09)
                ion.rotation_frequency_tor = ω
                any_ion_rot = true
            end
        end
        if any_ion_rot
            IMAS.ωtor2sonic!(cp1d)
        end
        @info "ActorPedestal: nn_predictor T_rot_ped = $(round(T_rot_krads; digits=3)) krad/s"
    end

    return actor
end

"""
    run_selected_pedestal_model(actor::ActorPedestal; density_factor::Float64, zeff_factor::Float64)

Runs selected pedestal model this prevents code duplication for using different par.model settings
"""
function run_selected_pedestal_model(actor::ActorPedestal; density_factor::Float64, zeff_factor::Float64)
    dd = actor.dd
    par = actor.par

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    if par.density_match == :ne_ped
        pedestal_density_tanh(dd, par; density_factor, zeff_factor, nn_prediction=actor.nn_prediction)
        finalize(step(actor.ped_actor))

    elseif par.density_match == :ne_line
        # NOTE: All pedestal actors take ne_ped as input
        # Here we convert the desirred pulse_schedule ne_line to ne_ped
        @assert par.ne_from in (:pulse_schedule, :nn_predictor) ":ne_line density_match requires ne_from ∈ (:pulse_schedule, :nn_predictor)"

        # run pedestal model on scaled density
        par.ne_from = :core_profiles
        pedestal_density_tanh(dd, par; density_factor=1.0, zeff_factor, nn_prediction=actor.nn_prediction)

        try
            # scale thermal densities to match desired line average (and temperatures accordingly, in case they matter)
            # we can do this because EPED and WPED only operate on temperature profiles
            nel_wanted = IMAS.ne_line(dd.pulse_schedule) * density_factor
            nel = IMAS.ne_line(eqt, cp1d)
            factor = nel_wanted / nel
            cp1d.electrons.density_thermal = cp1d.electrons.density_thermal * factor
            IMAS.unfreeze!(cp1d.electrons, :density)
            for ion in cp1d.ion
                if !ismissing(ion, :density_thermal)
                    ion.density_thermal = ion.density_thermal * factor
                    IMAS.unfreeze!(ion, :density)
                end
            end
            cp1d.electrons.temperature = cp1d.electrons.temperature / factor
            for ion in cp1d.ion
                if !ismissing(ion, :temperature)
                    ion.temperature = ion.temperature / factor
                end
            end

            # run the pedestal model
            finalize(step(actor.ped_actor))

        finally
            par.ne_from = :pulse_schedule
        end

    else
        error("act.ActorPedestal.density_match can be either one of [:ne_ped, :ne_line]")
    end

    # Apply NN-predicted Te / Ti / T_rot pedestal values on top of whatever
    # the underlying ped_actor produced. No-op unless par.ne_from == :nn_predictor.
    pedestal_nn_apply!(actor)

    return actor
end

function ti_te_ratio(cp1d, T_ratio_pedestal, rho_nml, rho_ped)
    if T_ratio_pedestal == 0.0
        # take ratio inside of the plasma core
        return IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average ./ cp1d.electrons.temperature)(rho_nml - (rho_ped - rho_nml))
    elseif T_ratio_pedestal <= 0.0
        return IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average ./ cp1d.electrons.temperature)(abs(T_ratio_pedestal))
    else
        return T_ratio_pedestal
    end
end

function _step(replay_actor::ActorReplay, actor::ActorPedestal, replay_dd::IMAS.dd)
    dd = actor.dd
    par = actor.par

    time0 = dd.global_time
    cp1d = dd.core_profiles.profiles_1d[time0]
    replay_cp1d = replay_dd.core_profiles.profiles_1d[time0]
    rho = cp1d.grid.rho_tor_norm

    # densities
    cp1d.electrons.density_thermal = IMAS.blend_core_edge(cp1d.electrons.density_thermal, replay_cp1d.electrons.density_thermal, rho, par.rho_nml, par.rho_ped; method=:shift)
    IMAS.unfreeze!(cp1d.electrons, :density)
    for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
        if !ismissing(ion, :density)
            ion.density = IMAS.blend_core_edge(ion.density, replay_ion.density, rho, par.rho_nml, par.rho_ped; method=:shift)
            IMAS.unfreeze!(ion, :density_thermal)
            if IMAS.hasdata(ion, :density_fast)
                ion.density_fast .= min.(ion.density_fast, ion.density)  # can't have more fast than total
            end
        end
    end

    # temperatures
    cp1d.electrons.temperature = IMAS.blend_core_edge(cp1d.electrons.temperature, replay_cp1d.electrons.temperature, rho, par.rho_nml, par.rho_ped)
    for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
        if !ismissing(ion, :temperature)
            ion.temperature = IMAS.blend_core_edge(ion.temperature, replay_ion.temperature, rho, par.rho_nml, par.rho_ped)
        end
    end

    # rotation (core from simulation, edge from replay)
    # NOTE: Cannot use blend_core_edge for rotation because it uses log internally,
    # and rotation can be negative. Use manual index-based shifting instead.
    i_nml = IMAS.argmin_abs(rho, par.rho_nml)
    ω_core = cp1d.rotation_frequency_tor_sonic
    ω_edge = replay_cp1d.rotation_frequency_tor_sonic
    ω_core[i_nml+1:end] = ω_edge[i_nml+1:end]
    ω_core[1:i_nml] = ω_core[1:i_nml] .- ω_core[i_nml] .+ ω_edge[i_nml]
    cp1d.rotation_frequency_tor_sonic = ω_core

    # Ion rotation: same blending
    for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
        if IMAS.hasdata(replay_ion, :rotation_frequency_tor)
            ω_ion_core = ion.rotation_frequency_tor
            ω_ion_edge = replay_ion.rotation_frequency_tor
            ω_ion_core[i_nml+1:end] = ω_ion_edge[i_nml+1:end]
            ω_ion_core[1:i_nml] = ω_ion_core[1:i_nml] .- ω_ion_core[i_nml] .+ ω_ion_edge[i_nml]
            ion.rotation_frequency_tor = ω_ion_core
        end
    end

    return replay_actor
end
