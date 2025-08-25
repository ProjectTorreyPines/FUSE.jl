#= ============= =#
#  ActorPedestal  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPedestal{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
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
    rotation_model::Switch{Symbol} = Switch{Symbol}([:linear, :replay, :none], "-", "Rotation edge model"; default=:none)
    #== L to H and H to L transition model ==#
    tau_t::Entry{T} = Entry{T}("s", "Edge temperature LH transition tanh evolution time (95% of full transition)")
    tau_n::Entry{T} = Entry{T}("s", "Edge density LH transition tanh evolution time (95% of full transition)")
    density_ratio_L_over_H::Entry{T} = Entry{T}("-", "n_Lmode / n_Hmode")
    zeff_ratio_L_over_H::Entry{T} = Entry{T}("-", "zeff_Lmode / zeff_Hmode")
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
end

"""
    ActorPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the pedestal boundary condition (height and width)
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
    actor = ActorPedestal(dd, par, act, noop, wped_actor, eped_actor, analytic_actor, noop, noop, Symbol[], -Inf, -Inf, -Inf, IMAS.core_profiles__profiles_1d{D}())
    actor.replay_actor = ActorReplay(dd, act.ActorReplay, actor)
    return actor
end

"""
    _step(actor::ActorPedestal)

Runs actors to evaluate profiles at the edge of the plasma
"""
function _step(actor::ActorPedestal{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    if !ismissing(par, :mode_transitions)
        causal_transition_time = IMAS.nearest_causal_time(sort!(collect(keys(par.mode_transitions))), dd.global_time).causal_time
        mode = par.mode_transitions[causal_transition_time]
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
            @assert par.ne_from == :pulse_schedule ":dynamic pedestal model requires `act.ActorPedestal.ne_from = :pulse_schedule`"
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
            time0 = dd.global_time
            rho = cp1d.grid.rho_tor_norm
            replay_cp1d = actor.replay_actor.replay_dd.core_profiles.profiles_1d[time0]
            cp1d.rotation_frequency_tor_sonic =
                IMAS.blend_core_edge(cp1d.rotation_frequency_tor_sonic, replay_cp1d.rotation_frequency_tor_sonic, rho, par.rho_nml, par.rho_ped; method=:shift)
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
function pedestal_density_tanh(dd::IMAS.dd, par::OverrideParameters{P,FUSEparameters__ActorPedestal{P}}; density_factor::Float64, zeff_factor::Float64) where {P<:Real}
    cp1d = dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm

    # Throughout FUSE, the "pedestal" values are defined at rho=0.9
    rho09 = 0.9

    # density pedestal width to match the existing temperature pedestal width
    w_ped = IMAS.pedestal_tanh_width_half_maximum(rho, cp1d.electrons.temperature)

    ne_old = copy(cp1d.electrons.density_thermal)
    ne_ped = IMAS.get_from(dd, Val(:ne_ped), par.ne_from, rho09) * density_factor
    cp1d.electrons.density_thermal[end] = ne_ped / 4.0
    ne = IMAS.blend_core_edge_Hmode(cp1d.electrons.density_thermal, rho, ne_ped, w_ped, par.rho_nml, par.rho_ped; method=:scale)
    cp1d.electrons.density_thermal = ne = IMAS.ped_height_at_09(rho, ne, ne_ped)
    ratio = ne ./ ne_old

    for ion in cp1d.ion
        if !ismissing(ion, :density_thermal)
            ion.density_thermal = ion.density_thermal .* ratio
            ni_ped = IMAS.interp1d(rho, ion.density_thermal).(rho09)
            ion.density_thermal[end] = ni_ped / 4.0
            ni = IMAS.blend_core_edge_Hmode(ion.density_thermal, rho, ni_ped, w_ped, par.rho_nml, par.rho_ped; method=:scale)
            ion.density_thermal = IMAS.ped_height_at_09(rho, ni, ni_ped)
        end
    end

    #NOTE: Zeff can change after a pedestal actor is run, even though actors like EPED and WPED only operate on the temperature profiles.
    # This is because in FUSE the calculation of Zeff is temperature dependent.
    zeff_ped = IMAS.get_from(dd, Val(:zeff_ped), par.zeff_from, rho09) * zeff_factor
    IMAS.scale_ion_densities_to_target_zeff!(cp1d, rho09, zeff_ped)

    return nothing
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
        pedestal_density_tanh(dd, par; density_factor, zeff_factor)
        finalize(step(actor.ped_actor))

    elseif par.density_match == :ne_line
        # NOTE: All pedestal actors take ne_ped as input
        # Here we convert the desirred pulse_schedule ne_line to ne_ped
        @assert par.ne_from == :pulse_schedule

        # run pedestal model on scaled density
        par.ne_from = :core_profiles
        pedestal_density_tanh(dd, par; density_factor=1.0, zeff_factor)

        try
            # scale thermal densities to match desired line average (and temperatures accordingly, in case they matter)
            # we can do this because EPED and WPED only operate on temperature profiles
            nel_wanted = IMAS.ne_line(dd.pulse_schedule) * density_factor
            nel = IMAS.ne_line(eqt, cp1d)
            factor = nel_wanted / nel
            cp1d.electrons.density_thermal = cp1d.electrons.density_thermal * factor
            for ion in cp1d.ion
                if !ismissing(ion, :density_thermal)
                    ion.density_thermal = ion.density_thermal * factor
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
    for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
        if !ismissing(ion, :density_thermal)
            ion.density_thermal = IMAS.blend_core_edge(ion.density_thermal, replay_ion.density_thermal, rho, par.rho_nml, par.rho_ped; method=:shift)
        end
    end

    # temperatures
    cp1d.electrons.temperature = IMAS.blend_core_edge(cp1d.electrons.temperature, replay_cp1d.electrons.temperature, rho, par.rho_nml, par.rho_ped)
    for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
        if !ismissing(ion, :temperature)
            ion.temperature = IMAS.blend_core_edge(ion.temperature, replay_ion.temperature, rho, par.rho_nml, par.rho_ped)
        end
    end

    # rotation
    cp1d.rotation_frequency_tor_sonic =
        IMAS.blend_core_edge(cp1d.rotation_frequency_tor_sonic, replay_cp1d.rotation_frequency_tor_sonic, rho, par.rho_nml, par.rho_ped; method=:shift)

    return replay_actor
end
