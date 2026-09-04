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
    nn_ped_quantities::Switch{Symbol} = Switch{Symbol}(
        [:ne_lh, :all], "-",
        "NN predictor outputs to apply: `:ne_lh` (ne_ped and L/H classification only, Te/Ti from EPED/WPED) or `:all` (ne_ped, te_ped, ti_ped, and rotation from NN)";
        default=:ne_lh)
    nn_hmode_enter::Entry{T} = Entry{T}("-", "fuse29 H-mode gate: enter H-mode when the predicted hmode probability rises above this (Schmitt trigger, see docs/GATING.md of the model)"; default=0.7)
    nn_hmode_exit::Entry{T} = Entry{T}("-", "fuse29 H-mode gate: leave H-mode when the predicted hmode probability falls below this (must be <= nn_hmode_enter)"; default=0.3)
    nn_warmup_time::Entry{Float64} = Entry{Float64}("s", "fuse29 recurrent state starts from zero (shot start). If the simulation starts mid-shot, replay the actuators from this time on the model's 50 ms grid before the first prediction so the state is warmed up (~1 s is enough). Missing: start from a zero state at the first call.")
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorPedestal{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.DD{D}
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
    nn_predictor::Union{Nothing,Fuse29NN}
    nn_prediction::Union{Nothing,NamedTuple}
end

"""
    ActorPedestal(dd::IMAS.DD, act::ParametersAllActors; kw...)

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
function ActorPedestal(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorPedestal(dd, act.ActorPedestal, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorPedestal(dd::IMAS.DD{D}, par::FUSEparameters__ActorPedestal{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
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

    # Run the fuse29 pedestal predictor once per time step. One call provides
    # ne, Te, Ti and T_rot at rho_tor = 0.85 plus the H-mode probability; the
    # recurrent state is kept on the model's own 50 ms grid by fuse29_predict!
    # (exactly one model tick per 50 ms of simulated time, however often the
    # actor is evaluated within a step).
    if par.ne_from == :nn_predictor
        fuse29_predict!(actor)
    end
    if !ismissing(par, :mode_transitions)
        causal_transition_time = IMAS.nearest_causal_time(sort!(collect(keys(par.mode_transitions))), dd.global_time).causal_time
        mode = par.mode_transitions[causal_transition_time]
    elseif par.ne_from == :nn_predictor && actor.nn_prediction !== nothing
        nnp = actor.nn_prediction
        @info "t = $(dd.global_time) s | fuse29: hmode_prob = $(round(nnp.hmode_prob; digits=3)) -> $(nnp.is_h_mode ? "H" : "L")-mode, ne(ρ=$(nnp.rho_ref)) = $(round(nnp.ne; digits=3))e19 m⁻³, ticks = $(nnp.n_ticks)"
        # Hysteretic gate decided in fuse29_advance! (enter > nn_hmode_enter, exit < nn_hmode_exit)
        mode = nnp.is_h_mode ? :H_mode : :L_mode
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
                # save initial L-mode ratio for use in future H→L back-transitions

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
    pedestal_density_tanh(dd::IMAS.DD, par::OverrideParameters{P,FUSEparameters__ActorPedestal{P}}; density_factor::Float64, zeff_factor::Float64) where {P<:Real}

The edge density must be defined independently of the pedestal model

The EPED and WPED models only operate on the temperature profiles
"""
function pedestal_density_tanh(dd::IMAS.DD, par::OverrideParameters{P,FUSEparameters__ActorPedestal{P}};
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
        @assert nn_prediction !== nothing "ne_from=:nn_predictor requires nn_prediction (call fuse29_predict! first)"
        # fuse29 `ne` is the electron density AT rho_tor = nn_prediction.rho_ref (0.85),
        # valid in every regime — pin the profile there rather than at rho=0.9.
        ne_ped = Float64(nn_prediction.ne) * 1e19 * density_factor
        rho_pin = nn_prediction.rho_ref
        @info "ActorPedestal: nn_predictor ne(ρ=$(rho_pin)) = $(round(nn_prediction.ne * density_factor; digits=3)) x 10^19 m^-3"
    else
        ne_ped = IMAS.get_from(dd, Val(:ne_ped), par.ne_from, rho09) * density_factor
        rho_pin = rho09
    end
    cp1d.electrons.density_thermal[end] = ne_ped / 4.0
    ne = IMAS.blend_core_edge_Hmode(cp1d.electrons.density_thermal, rho, ne_ped, w_ped, par.rho_nml, par.rho_ped; method=:scale)
    cp1d.electrons.density_thermal = ne = _scale_profile_at(rho, ne, rho_pin, ne_ped)
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
    _scale_profile_at(rho, profile, rho_at, value)

Scale `profile` so that it equals `value` at `rho_at` (generalises
`IMAS.ped_height_at_09`, which pins at rho = 0.9).
"""
function _scale_profile_at(rho::AbstractVector{<:Real}, profile::AbstractVector{<:Real}, rho_at::Real, value::Real)
    return profile ./ IMAS.interp1d(rho, profile).(rho_at) .* value
end

"""
    build_fuse29_actuators(nn::Fuse29NN, dd::IMAS.DD; source::Symbol=:zmq, time::Float64=dd.global_time) -> (u::Vector{Float32}, live::Vector{Bool})

Assemble the 29-channel fuse29 actuator vector, in **raw physical units** and
in `FUSE29_ACTUATOR_NAMES` order, for the 50 ms interval ending at `time`.
Channels with no source are left at the training-split **median** (the model
authors' recommended stand-in; the mean vector is jointly unphysical and zeros
put the coil channels several sigma out) and their `live` flag stays `false`.

`source = :zmq` reads the live GSLite signals stored in `dd._aux` by
`ActorZMQ.receive!` (latest causal sample at `time`); `source = :dd` reads
standard IMAS fields (no ZMQ required).

| channel | units | `:zmq` | `:dd` |
|---|---|---|---|
| `pinj` | MW | `dd._aux[:zmq_Pnbi]` (W) / 1e6 | Σ `dd.pulse_schedule.nbi.unit[].power.reference` (W) / 1e6 |
| `tinj` | N·m | Σ `dd.core_sources` NBI (`identifier.index == 2`) `global_quantities[time].torque_tor` | same |
| `ech_total` | MW † | `dd._aux[:zmq_Pech]` (W) / 1e6 | Σ `dd.ec_launchers.beam[].power_launched` (W) / 1e6 |
| `f1a..f9b` | A | `dd._aux[:zmq_I_coil][7..24]` (PCF1A..PCF9B) | `dd.pf_active.coil` named `F1A..F9B`, current / `turns_with_sign` |
| `ecoila`, `ecoilb` | A | `dd._aux[:zmq_I_coil][1]`, `[4]` (PCECOILA, PCECOILB) | `dd.pf_active.coil` named `ECOILA`, `ECOILB` |
| `gasa_cal..gase_cal` | Torr·L/s | `dd._aux[:zmq_gas[a-e]_cal]` | (no dd source → median) |
| `bt` | T (signed) | `dd.equilibrium.vacuum_toroidal_field.b0` at `time` | same |

`I_coil` is the 24-element PCS coil-current vector (PCECOILA, PCE89DN, PCE567UP,
PCECOILB, PCE89UP, PCE567DN, PCF1A..PCF9A, PCF1B..PCF9B); the internal C-coil
segments at 2,3,5,6 are not model inputs.

† `ech_total` carries a scale factor inherited from the upstream DIII-D
pipeline: its training p99 is ~30 "MW", above the installed ECH power, and the
median is 0. A nominal-MW value sits inside the training distribution, so pass
MW and let `fuse29_check_units` judge; treat a mismatch on this channel as
expected rather than as evidence the wiring is wrong.

The predecessor's `pohm`, `ip` and `ipspr15v` inputs are gone: they are plasma
responses, not commands, and fuse29 deliberately does not take them.
"""
function build_fuse29_actuators(nn::Fuse29NN, dd::IMAS.DD; source::Symbol=:zmq, time::Float64=dd.global_time)
    u = fuse29_median_actuators(nn)
    live = falses(FUSE29_N_ACTUATORS)
    aux = getfield(dd, :_aux)

    function _set!(name::AbstractString, value)
        value === nothing && return
        isfinite(value) || return
        idx = findfirst(==(name), nn.actuator_names)
        idx === nothing && error("build_fuse29_actuators: channel \"$name\" is not in the bundle contract")
        u[idx] = Float32(value)
        live[idx] = true
        return
    end

    # Latest causal sample (t_i <= time); fall back to the first sample if all
    # entries are in the future (e.g. just after a time-rewind).
    function _aux_value_at(key::Symbol)
        haskey(aux, key) || return nothing
        rec = aux[key]
        (hasproperty(rec, :times) && hasproperty(rec, :values)) || return nothing
        isempty(rec.times) && return nothing
        idx = findlast(τ -> τ <= time + 1e-9, rec.times)
        idx === nothing && (idx = 1)
        return rec.values[idx]
    end

    function _interp_at(t::AbstractVector, y::AbstractVector, scheme::Symbol)
        (isempty(t) || length(t) != length(y)) && return nothing
        length(t) == 1 && return y[1]
        return IMAS.interp1d(t, y, scheme)(time)
    end

    if source == :dd
        # pinj — total NBI power from pulse_schedule (W -> MW)
        if !isempty(dd.pulse_schedule.nbi.unit) && !isempty(dd.pulse_schedule.nbi.time)
            Pnbi = 0.0
            found = false
            for unit in dd.pulse_schedule.nbi.unit
                if !ismissing(unit.power, :reference) && !isempty(unit.power.reference)
                    v = _interp_at(dd.pulse_schedule.nbi.time, unit.power.reference, :constant)
                    v === nothing && continue
                    Pnbi += v
                    found = true
                end
            end
            found && _set!("pinj", Pnbi / 1e6)
        end

        # ech_total — total EC power from ec_launchers (W -> "MW", see docstring)
        if !isempty(dd.ec_launchers.beam)
            Pech = 0.0
            found = false
            for beam in dd.ec_launchers.beam
                if !ismissing(beam, :power_launched) && !isempty(beam.power_launched.time)
                    v = _interp_at(beam.power_launched.time, beam.power_launched.data, :constant)
                    v === nothing && continue
                    Pech += v
                    found = true
                end
            end
            found && _set!("ech_total", Pech / 1e6)
        end

        # coil currents from dd.pf_active — DIII-D names as loaded from the OMFIT/D3D machine mapping
        coil_name_map = Dict(
            "ECOILA" => "ecoila", "ECOILB" => "ecoilb",
            "F1A" => "f1a", "F2A" => "f2a", "F3A" => "f3a", "F4A" => "f4a", "F5A" => "f5a",
            "F6A" => "f6a", "F7A" => "f7a", "F8A" => "f8a", "F9A" => "f9a",
            "F1B" => "f1b", "F2B" => "f2b", "F3B" => "f3b", "F4B" => "f4b", "F5B" => "f5b",
            "F6B" => "f6b", "F7B" => "f7b", "F8B" => "f8b", "F9B" => "f9b"
        )
        for coil in dd.pf_active.coil
            ch_name = get(coil_name_map, uppercase(strip(coil.name)), nothing)
            ch_name === nothing && continue
            if !ismissing(coil.current, :data) && !isempty(coil.current.data)
                turns = isempty(coil.element) ? 1.0 : coil.element[1].turns_with_sign
                v = _interp_at(coil.current.time, coil.current.data, :linear)
                v === nothing && continue
                _set!(ch_name, v / turns)
            end
        end
        # gasa..gase_cal — no dd source; stay at the training median

    elseif source == :zmq
        let v = _aux_value_at(:zmq_Pnbi); v === nothing || _set!("pinj", v / 1e6); end
        let v = _aux_value_at(:zmq_Pech); v === nothing || _set!("ech_total", v / 1e6); end
        for (k, name) in zip(
                (:zmq_gasa_cal, :zmq_gasb_cal, :zmq_gasc_cal, :zmq_gasd_cal, :zmq_gase_cal),
                ("gasa_cal",    "gasb_cal",    "gasc_cal",    "gasd_cal",    "gase_cal"))
            v = _aux_value_at(k)
            v === nothing || _set!(name, v)
        end
        let v = _aux_value_at(:zmq_I_coil)
            if v !== nothing
                length(v) >= 1 && _set!("ecoila", v[1])
                length(v) >= 4 && _set!("ecoilb", v[4])
                f_names = ("f1a", "f2a", "f3a", "f4a", "f5a", "f6a", "f7a", "f8a", "f9a",
                           "f1b", "f2b", "f3b", "f4b", "f5b", "f6b", "f7b", "f8b", "f9b")
                for (k, name) in enumerate(f_names)
                    length(v) >= 6 + k || break
                    _set!(name, v[6+k])
                end
            end
        end

    else
        error("build_fuse29_actuators: source must be :zmq or :dd, got $(repr(source))")
    end

    # tinj — NBI torque from core_sources for both paths (NBI identifier index = 2)
    tinj = 0.0
    found_tinj = false
    for src in dd.core_sources.source
        if src.identifier.index == 2 && !isempty(src.global_quantities)
            gq = try
                src.global_quantities[time]
            catch
                src.global_quantities[end]
            end
            if !ismissing(gq, :torque_tor)
                tinj += gq.torque_tor
                found_tinj = true
            end
        end
    end
    found_tinj && _set!("tinj", tinj)

    # bt — signed vacuum toroidal field at `time` for both paths
    # (b0's time coordinate is dd.equilibrium.time, resolved by get_time_array)
    let vtf = dd.equilibrium.vacuum_toroidal_field
        if !ismissing(vtf, :b0) && !isempty(vtf.b0)
            v = try
                IMAS.get_time_array(vtf, :b0, time, :constant)
            catch
                time <= dd.equilibrium.time[1] ? vtf.b0[1] : vtf.b0[end]
            end
            _set!("bt", v)
        end
    end

    return u, live
end

"""
    fuse29_predict!(actor::ActorPedestal)

Advance the fuse29 pedestal predictor to `dd.global_time` and store the result
in `actor.nn_prediction`. Loads the model lazily (cached per process) and keeps
the per-shot recurrent state in `dd._aux[:fuse29]` (a [`Fuse29Tracker`](@ref)),
so the model sees exactly one tick per 50 ms of simulated time regardless of
how often the actor is evaluated: catch-up ticks when the host step is longer
than 50 ms, zero-order hold when it is shorter, re-step from the last committed
state on a repeated/iterated step, and a fresh zero state (with a warning) when
time is rewound. `par.nn_warmup_time` replays the actuators from that time on
the first call so a mid-shot start does not begin from a cold state (the model
needs ~1 s to settle from zeros).

`nn_prediction` fields: `ne` (1e19 m⁻³), `te_ped`, `ti_ped` (keV), `t_rot_ped`
(krad/s) — profile values at `rho_ref` = 0.85, valid in every regime;
`neped_prmtan`, `teped_prmtan`, `rho_sym`, `ne_top_loc` — tanh-fit pedestal
height/location, `NaN` unless `is_h_mode`; `hmode_prob`, `is_h_mode` (Schmitt
trigger with `nn_hmode_enter`/`nn_hmode_exit`), `raw`, `rho_ref`, `time`, `n_ticks`.
"""
function fuse29_predict!(actor::ActorPedestal)
    dd = actor.dd
    par = actor.par

    if actor.nn_predictor === nothing
        actor.nn_predictor = load_pedestal_nn()
        @info "ActorPedestal: loaded fuse29 pedestal predictor $(actor.nn_predictor)"
    end
    nn = actor.nn_predictor
    aux = getfield(dd, :_aux)
    t_now = dd.global_time
    eps = 1e-3 * nn.period_s

    tr = get(aux, :fuse29, nothing)
    if !(tr isa Fuse29Tracker) || t_now < tr.committed_time - eps
        if tr isa Fuse29Tracker
            @warn "ActorPedestal: fuse29 time rewound from $(tr.committed_time) s to $t_now s; restarting the recurrent state from zero"
        end
        t_start = t_now - nn.period_s
        if !ismissing(par, :nn_warmup_time)
            t_start = min(par.nn_warmup_time, t_start)
        end
        tr = Fuse29Tracker(nn, t_start)
        aux[:fuse29] = tr
    end

    function actuators_at(t::Float64)
        u, live = build_fuse29_actuators(nn, dd; source=par.fpe_source, time=t)
        if !tr.units_checked
            tr.units_checked = true
            for w in fuse29_check_units(nn, u)
                @warn "ActorPedestal: fuse29 unit check — $w"
            end
            missing_ch = [String(name) for (name, l) in zip(nn.actuator_names, live) if !l]
            if !isempty(missing_ch)
                @warn "ActorPedestal: fuse29 channels with no `$(par.fpe_source)` source, held at the training median: $(join(missing_ch, ", "))"
            end
        end
        exc = fuse29_in_distribution(nn, u)
        isempty(exc) || @debug "ActorPedestal: fuse29 actuators outside training p1..p99 at t=$t s" excursions = exc
        return u
    end

    pred, is_h, n = fuse29_advance!(tr, nn, t_now, actuators_at;
                                    hmode_enter=Float64(par.nn_hmode_enter), hmode_exit=Float64(par.nn_hmode_exit))
    actor.nn_prediction = fuse29_prediction(nn, pred, is_h; time=t_now, n_ticks=n)
    return actor
end

"""
    pedestal_nn_apply!(actor::ActorPedestal)

When `par.ne_from == :nn_predictor` and `par.nn_ped_quantities == :all`, blend
the fuse29-predicted edge temperatures and rotation into `cp1d` on top of
whatever the underlying pedestal actor produced, so downstream consumers (and
`_finalize`'s `dd.summary.local.pedestal` writes) reflect them.

- `nn_prediction.te_ped`    (keV at rho_ref) -> `cp1d.electrons.temperature`
- `nn_prediction.ti_ped`    (keV at rho_ref) -> every `cp1d.ion[*].temperature`
- `nn_prediction.t_rot_ped` (krad/s at rho_ref) -> every `cp1d.ion[*].rotation_frequency_tor`,
   ONLY when `par.rotation_model == :nn_pedestal` (explicit opt-in). The default
   `:none` leaves rotation untouched.

These heads are profile values at `rho_ref` (0.85) and are valid in every
regime: in H-mode the profile is given a tanh pedestal and pinned to the NN
value at `rho_ref`; in L-mode the edge is scaled WPED-style
(`IMAS.blend_core_edge_Lmode`) to hit the NN value at `rho_ref`, without
manufacturing a pedestal.
"""
function pedestal_nn_apply!(actor::ActorPedestal)
    par = actor.par
    nn = actor.nn_prediction
    (par.ne_from == :nn_predictor && nn !== nothing) || return actor

    # When :ne_lh, Te/Ti/rotation come from EPED/WPED; NN used only for ne and L/H.
    if par.nn_ped_quantities == :ne_lh
        return actor
    end

    is_h_mode = !isempty(actor.state) && actor.state[end] == :H_mode
    cp1d = actor.dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm
    rho_ref = nn.rho_ref
    w_ped = IMAS.pedestal_tanh_width_half_maximum(rho, cp1d.electrons.temperature)

    # Snapshot Ti/Te ratio *before* we mutate Te so the Ti separatrix boundary
    # stays consistent with the post-ped_actor profiles.
    Ti_over_Te = ti_te_ratio(cp1d, par.T_ratio_pedestal, par.rho_nml, par.rho_ped)

    function _apply_edge(profile::AbstractVector, value_eV::Float64, sep::Float64)
        prof = copy(profile)
        prof[end] = sep
        if is_h_mode
            prof = IMAS.blend_core_edge_Hmode(prof, rho, value_eV, w_ped, par.rho_nml, par.rho_ped; method=:scale)
            return _scale_profile_at(rho, prof, rho_ref, value_eV)
        else
            return IMAS.blend_core_edge_Lmode(prof, rho, value_eV, rho_ref)
        end
    end

    Te_ped_keV = nn.te_ped
    if isfinite(Te_ped_keV) && Te_ped_keV > 0
        cp1d.electrons.temperature = _apply_edge(cp1d.electrons.temperature, Te_ped_keV * 1e3, Float64(par.Te_sep))
    end

    Ti_ped_keV = nn.ti_ped
    if isfinite(Ti_ped_keV) && Ti_ped_keV > 0
        Ti_sep = Float64(par.Te_sep) * Ti_over_Te
        for ion in cp1d.ion
            if !ismissing(ion, :temperature)
                ion.temperature = _apply_edge(ion.temperature, Ti_ped_keV * 1e3, Ti_sep)
            end
        end
    end

    # rotation_frequency_tor can be negative, so we cannot use blend_core_edge_Hmode
    # (which uses log internally — see the :replay rotation branch in _step).
    # Use a simple scale-to-target at rho_ref, falling back to a constant offset
    # if the existing trace passes through zero there. We also refresh
    # cp1d.rotation_frequency_tor_sonic so FINN's
    # `profile_from_rotation_shear_transport` sees the NN pedestal value as BC.
    T_rot_krads = nn.t_rot_ped
    if isfinite(T_rot_krads) && par.rotation_model == :nn_pedestal
        ω_ped = T_rot_krads * 1e3
        any_ion_rot = false
        for ion in cp1d.ion
            if !ismissing(ion, :rotation_frequency_tor)
                ω = copy(ion.rotation_frequency_tor)
                ω_at_ref = IMAS.interp1d(rho, ω).(rho_ref)
                ω = iszero(ω_at_ref) ? ω .+ ω_ped : ω .* (ω_ped / ω_at_ref)
                ion.rotation_frequency_tor = ω
                any_ion_rot = true
            end
        end
        if any_ion_rot
            IMAS.ωtor2sonic!(cp1d)
        end
        @info "ActorPedestal: nn_predictor T_rot(ρ=$(rho_ref)) = $(round(T_rot_krads; digits=3)) krad/s"
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

        # save original ne_from so we can restore it (not hardcode :pulse_schedule)
        original_ne_from = par.ne_from

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
            par.ne_from = original_ne_from
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

function _step(replay_actor::ActorReplay, actor::ActorPedestal, replay_dd::IMAS.DD)
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
