import EPEDNN

#= ============= =#
#  ActorPedestal  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorPedestal{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
    density_match::Switch{Symbol} = Switch{Symbol}([:ne_line, :ne_ped], "-", "Matching density based on ne_ped or line averaged density"; default=:ne_ped)
    model::Switch{Symbol} = Switch{Symbol}([:EPED, :WPED, :dynamic, :none], "-", "Pedestal model to use"; default=:EPED)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    ne_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_ped_from::Switch{Symbol} = switch_get_from(:zeff_ped)

    #== L to H and H to L transition model ==#
    tau_t::Entry{T} = Entry{T}("[s]", "pedestal temperature LH transition tanh evolution time (95% of full transition)")
    tau_n::Entry{T} = Entry{T}("[s]", "pedestal density LH transition tanh evolution time (95% of full transition)")
    density_ratio_L_over_H::Entry{T} = Entry{T}("[-]", "n_Lmode / n_Hmode")
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
    eped_actor = ActorEPED(dd, act.ActorEPED; ne_ped_from=par.ne_from, par.zeff_ped_from, par.βn_from, par.ip_from, par.rho_nml, par.rho_ped)
    wped_actor = ActorWPED(dd, act.ActorWPED; ne_ped_from=par.ne_from, par.zeff_ped_from, par.rho_ped, par.do_plot)
    none_actor = ActorNoOperation(dd, act.ActorNoOperation)
    return ActorPedestal(dd, par, act, none_actor, none_actor, eped_actor, wped_actor, Symbol[], -Inf, -Inf, IMAS.core_profiles__profiles_1d())
end

"""
    _step(actor::ActorPedestal)

Runs actors to evaluate profiles at the edge of the plasma
"""
function _step(actor::ActorPedestal{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    if par.model == :none
        actor.ped_actor = actor.none_actor
    elseif par.model == :EPED
        actor.ped_actor = actor.eped_actor
        run_selected_pedstal_model(actor)
    elseif par.model == :WPED
        actor.ped_actor = actor.wped_actor
        run_selected_pedstal_model(actor)
    elseif par.model == :dynamic
        @assert par.ne_from == :pulse_schedule "Dyanmic requires ne_from pulse schedule"

        if IMAS.satisfies_h_mode_conditions(dd; threshold_multiple=2.0)
            push!(actor.state,:H_mode)
        elseif !IMAS.satisfies_h_mode_conditions(dd; threshold_multiple=1.0)
            push!(actor.state,:L_mode)
        else
            push!(actor.state,actor.state[end])
        end

        if length(actor.state) >= 2 && actor.state[end-1:end] == [:L_mode, :H_mode] 
            # L to H
            actor.t_lh = dd.global_time
        elseif length(actor.state) >= 2 && actor.state[end-1:end] == [:H_mode, :L_mode] 
            # H to L
            actor.t_hl = dd.global_time

            # The L to H and H to L model are triggered when hmode conditions are triggered so therefore pulse_schedule needs to be updated based on this model
            if par.density_match == :ne_ped
                @ddtime(dd.pulse_schedule.density_control.n_e_pedestal.reference =  @ddtime(dd.summary.local.pedestal.n_e.value))
            elseif par.density_match == :ne_line
                @ddtime(dd.pulse_schedule.density_control.n_e_line.reference =  IMAS.geometric_midplane_line_averaged_density(eqt, cp1d))
            end
        end

        if actor.state[end] == :L_mode
            α_t = LH_tanh_response(par.tau_t, dd.global_time,actor.t_hl;time_steps=100)
            α_n = LH_tanh_response(par.tau_n, dd.global_time,actor.t_hl;time_steps=100)

            actor.ped_actor = actor.wped_actor
            actor.ped_actor.par.density_factor = 1. - α_n * (1 - par.density_ratio_L_over_H)

            run_selected_pedstal_model(actor)

            Te_ongoing = (1 .- α_t) .* actor.cp1d_transition.electrons.temperature .+  α_t .* dd.core_profiles.profiles_1d[].electrons.temperature
            Ti_ongoing = (1 .- α_t) .* actor.cp1d_transition.ion[1].temperature .+  α_t .* dd.core_profiles.profiles_1d[].ion[1].temperature

            cp1d.electrons.temperature = Te_ongoing
            for ion in cp1d.ion
                ion.temperature = Ti_ongoing
            end
        else
            # H mode
            α_t = LH_tanh_response(par.tau_t,dd.global_time,actor.t_lh;time_steps=100)
            α_n = LH_tanh_response(par.tau_n,dd.global_time, actor.t_lh;time_steps=100)
            actor.ped_actor = actor.eped_actor
           
            actor.ped_actor.par.density_factor =  1 + (1 / par.density_ratio_L_over_H - 1) * α_n
            run_selected_pedstal_model(actor)

            Te_ongoing = (1 .- α_t) .* actor.cp1d_transition.electrons.temperature .+  α_t .* dd.core_profiles.profiles_1d[].electrons.temperature
            Ti_ongoing = (1 .- α_t) .* actor.cp1d_transition.ion[1].temperature .+  α_t .* dd.core_profiles.profiles_1d[].ion[1].temperature

            cp1d.electrons.temperature = Te_ongoing
            for ion in cp1d.ion
                ion.temperature = Ti_ongoing
            end
        end
    end

    return actor
end

"""
    LH_tanh_response(τ::Float64,t_now::Float64, t_LH::Float64; time_steps::Int=100)

Returns a parameter that follows a tanh like response where τ represent the value of 0.95 @ τ time starting from t_LH
"""
function LH_tanh_response(τ::Float64,t_now::Float64, t_LH::Float64; time_steps::Int=100)
    if isinf(t_LH)
        return 0.0
    end
    t, α = LH_tanh_response(τ, t_LH; time_steps)
    return maximum([minimum([IMAS.interp1d(t,α).(t_now),1.0]),0.0])
end

function LH_tanh_response(τ::Float64,t_LH::Float64; time_steps::Int=100)
    t = (0:τ/time_steps:τ) .+ t_LH
    α = tanh.((2pi.*(t .- t_LH .- τ / 4.0)) ./ τ)
    α = ((α .- α[1]) ./ (α[end]-α[1]))
    return t, α
end

"""
    run_selected_pedstal_model(actor::ActorPedestal)

Runs selected pedestal model this prevents code duplication for using different par.model settings
"""
function run_selected_pedstal_model(actor::ActorPedestal)
    par = actor.par
    if par.density_match == :ne_ped
        finalize(step(actor.ped_actor))

    elseif par.density_match == :ne_line
        # NOTE: All pedestal actors take ne_ped as input
        # Here we convert the desirred pulse_schedule ne_line to ne_ped
        @assert par.ne_from == :pulse_schedule

        # freeze pressures since they are input to the pedestal models
        IMAS.refreeze!(cp1d, :pressure_thermal)
        IMAS.refreeze!(cp1d, :pressure)

        # run pedestal model on scaled density
        if par.model != :none
            actor.ped_actor.par.ne_ped_from = :core_profiles
        end

        # first run the pedestal model on the density as is
        _finalize(_step(actor.ped_actor))

        # scale thermal densities to match desired line average
        ne_line = IMAS.geometric_midplane_line_averaged_density(eqt, cp1d)
        ne_line_wanted = IMAS.ne_line(dd.pulse_schedule)
        factor = ne_line_wanted / ne_line
        cp1d.electrons.density_thermal = cp1d.electrons.density_thermal * factor
        for ion in cp1d.ion
            ion.density_thermal = ion.density_thermal * factor
        end

        finalize(step(actor.ped_actor))
        if par.model != :none
            actor.ped_actor.par.ne_ped_from = :pulse_schedule
        end

        # turn pressures back into expressions
        IMAS.empty!(cp1d, :pressure_thermal)
        IMAS.empty!(cp1d, :pressure)
    end
end