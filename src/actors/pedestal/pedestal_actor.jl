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
    #== actor parameters==#
    density_match::Switch{Symbol} = Switch{Symbol}([:ne_line, :ne_ped], "-", "Matching density based on ne_ped or line averaged density"; default=:ne_ped)
    model::Switch{Symbol} = Switch{Symbol}([:EPED, :WPED, :dynamic, :none], "-", "Pedestal model to use"; default=:EPED)
    #== L to H and H to L transition model ==#
    tau_t::Entry{T} = Entry{T}("s", "pedestal temperature LH transition tanh evolution time (95% of full transition)")
    tau_n::Entry{T} = Entry{T}("s", "pedestal density LH transition tanh evolution time (95% of full transition)")
    density_ratio_L_over_H::Entry{T} = Entry{T}("-", "n_Lmode / n_Hmode")
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorPedestal{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPedestal{P}
    act::ParametersAllActors{P}
    ped_actor::Union{ActorNoOperation{D,P},ActorEPED{D,P},ActorWPED{D,P}}
    none_actor::ActorNoOperation{D,P}
    eped_actor::ActorEPED{D,P}
    wped_actor::ActorWPED{D,P}
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

function ActorPedestal(dd::IMAS.dd, par::FUSEparameters__ActorPedestal, act::ParametersAllActors; kw...)
    logging_actor_init(ActorPedestal)
    par = par(kw...)
    eped_actor = ActorEPED(dd, act.ActorEPED; par.rho_nml, par.rho_ped, par.T_ratio_pedestal, par.Te_sep, par.ip_from, par.βn_from, par.ne_from, par.zeff_from)
    wped_actor = ActorWPED(dd, act.ActorWPED; par.rho_nml, par.rho_ped, par.T_ratio_pedestal, par.Te_sep, par.ip_from, par.βn_from, par.ne_from, par.zeff_from)
    none_actor = ActorNoOperation(dd, act.ActorNoOperation)
    return ActorPedestal(dd, par, act, none_actor, none_actor, eped_actor, wped_actor, Symbol[], -Inf, -Inf, -Inf, IMAS.core_profiles__profiles_1d())
end

"""
    _step(actor::ActorPedestal)

Runs actors to evaluate profiles at the edge of the plasma
"""
function _step(actor::ActorPedestal{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    debug = false

    if par.model == :none
        actor.ped_actor = actor.none_actor
        return actor
    end

    # Make all densities conform to EPED tanh form
    # the EPED and WPED models only operate on the temperature profiles
    rho09 = 0.9 # FUSE defines "pedestal" as rho=0.9
    w_ped = 1.0 - rho09
    rho = cp1d.grid.rho_tor_norm
    ne_ped = IMAS.interp1d(rho, cp1d.electrons.density_thermal).(rho09)
    cp1d.electrons.density_thermal[end] = ne_ped / 4.0
    cp1d.electrons.density_thermal = IMAS.blend_core_edge_Hmode(cp1d.electrons.density_thermal, rho, ne_ped, w_ped, par.rho_nml, par.rho_ped)
    for ion in cp1d.ion
        if !ismissing(ion, :density_thermal)
            ni_ped = IMAS.interp1d(rho, ion.density_thermal).(rho09)
            ion.density_thermal[end] = ni_ped / 4.0
            ion.density_thermal = IMAS.blend_core_edge_Hmode(ion.density_thermal, rho, ni_ped, w_ped, par.rho_nml, par.rho_ped)
        end
    end

    if par.model == :EPED
        actor.ped_actor = actor.eped_actor
        run_selected_pedstal_model(actor)

    elseif par.model == :WPED
        actor.ped_actor = actor.wped_actor
        run_selected_pedstal_model(actor)

    elseif par.model == :dynamic
        @assert par.ne_from == :pulse_schedule ":dynamic pedestal model requires `act.ActorPedestal.ne_from = :pulse_schedule`"
        @assert actor.previous_time < dd.global_time "subsequent calls to :dynamic pedestal model require dd.global_time advance"

        if IMAS.satisfies_h_mode_conditions(dd; threshold_multiplier=1.5)
            push!(actor.state, :H_mode)
        elseif !IMAS.satisfies_h_mode_conditions(dd; threshold_multiplier=1.0)
            push!(actor.state, :L_mode)
        else
            push!(actor.state, actor.state[end])
        end

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

        if actor.state[end] == :L_mode
            # L-mode
            α_t = LH_tanh_response(par.tau_t, actor.t_hl, dd.global_time) # from 0 -> 1
            α_n = LH_tanh_response(par.tau_n, actor.t_hl, dd.global_time) # from 0 -> 1

            actor.ped_actor = actor.wped_actor
            actor.ped_actor.par.density_factor = 1.0 * (1 - α_n) + par.density_ratio_L_over_H * α_n

            run_selected_pedstal_model(actor)

            Te_ongoing = (1 .- α_t) .* actor.cp1d_transition.electrons.temperature .+ α_t .* dd.core_profiles.profiles_1d[].electrons.temperature
            Ti_ongoing = (1 .- α_t) .* actor.cp1d_transition.ion[1].temperature .+ α_t .* dd.core_profiles.profiles_1d[].ion[1].temperature

            cp1d.electrons.temperature = Te_ongoing
            for ion in cp1d.ion
                ion.temperature = Ti_ongoing
            end

        else
            # H-mode
            α_t = LH_tanh_response(par.tau_t, actor.t_lh, dd.global_time) # from 0 -> 1
            α_n = LH_tanh_response(par.tau_n, actor.t_lh, dd.global_time) # from 0 -> 1

            actor.ped_actor = actor.eped_actor
            actor.ped_actor.par.density_factor = par.density_ratio_L_over_H * (1 - α_n) + 1.0 * α_n

            run_selected_pedstal_model(actor)

            Te_ongoing = (1 .- α_t) .* actor.cp1d_transition.electrons.temperature .+ α_t .* dd.core_profiles.profiles_1d[].electrons.temperature
            Ti_ongoing = (1 .- α_t) .* actor.cp1d_transition.ion[1].temperature .+ α_t .* dd.core_profiles.profiles_1d[].ion[1].temperature

            cp1d.electrons.temperature = Te_ongoing
            for ion in cp1d.ion
                ion.temperature = Ti_ongoing
            end
        end

        if debug
            println()
            @show dd.global_time
            println(actor.state[end])
            @show actor.t_lh
            @show actor.t_hl
            @show α_t
            @show α_n
        end

        actor.previous_time = dd.global_time

    end

    return actor
end

"""
    LH_tanh_response(τ::Float64, t_LH::Float64, t_now::Float64)

Returns a parameter that follows a tanh like response where τ represent the value of 0.95 @ τ time starting from t_LH
"""
function LH_tanh_response(τ::Float64, t_LH::Float64, t_now::Float64)
    if t_LH == -Inf
        return 1.0
    elseif t_now <= t_LH
        return 0.0
    end
    α = tanh.((2pi .* (t_now .- t_LH .- τ / 4.0)) ./ τ) / 2.0 + 0.5
    α0 = tanh.((2pi .* (.-τ / 4.0)) ./ τ) / 2.0 + 0.5
    α = (α .- α0) ./ (1 - α0)
    return α
end

"""
    run_selected_pedstal_model(actor::ActorPedestal)

Runs selected pedestal model this prevents code duplication for using different par.model settings
"""
function run_selected_pedstal_model(actor::ActorPedestal)
    dd = actor.dd
    par = actor.par

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    if par.density_match == :ne_ped
        finalize(step(actor.ped_actor))

    elseif par.density_match == :ne_line
        # NOTE: All pedestal actors take ne_ped as input
        # Here we convert the desirred pulse_schedule ne_line to ne_ped
        @assert par.ne_from == :pulse_schedule

        # run pedestal model on scaled density
        if par.model != :none
            actor.ped_actor.par.ne_from = :core_profiles
        end

        # scale thermal densities to match desired line average (and temperatures accordingly, in case they matter)
        # we can do this because EPED and WPED only operate on temperature profiles
        ne_line_wanted = IMAS.ne_line(dd.pulse_schedule) * actor.ped_actor.par.density_factor
        ne_line = IMAS.geometric_midplane_line_averaged_density(eqt, cp1d)
        factor = ne_line_wanted / ne_line
        cp1d.electrons.density_thermal = cp1d.electrons.density_thermal * factor
        for ion in cp1d.ion
            ion.density_thermal = ion.density_thermal * factor
        end
        cp1d.electrons.temperature = cp1d.electrons.temperature / factor
        for ion in cp1d.ion
            ion.temperature = ion.temperature / factor
        end

        # run the pedestal model
        finalize(step(actor.ped_actor))

        if par.model != :none
            actor.ped_actor.par.ne_from = :pulse_schedule
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