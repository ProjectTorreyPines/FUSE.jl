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

    # tau_t
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
        @show "dynamic"
        if IMAS.satisfies_h_mode_conditions(dd)
            push!(actor.state,:H_mode)
        else
            push!(actor.state,:L_mode)
        end
        neped = IMAS.get_from(dd, Val{:ne_ped}, actor.par.ne_from, nothing)
        @show neped
        # 2 scenarios
        @show length(actor.state)
        if length(actor.state) >= 2 && actor.state[end-1:end] == [:L_mode, :H_mode] 
            # L to H
            actor.t_lh = dd.global_time
        elseif length(actor.state) >= 2 && actor.state[end-1:end] == [:H_mode, :L_mode] 
            # H to L
            actor.t_hl = dd.global_time
        end

        # actor.cp1d_transition.electrons.temperature
        # actor.cp1d_transition.electrons.density_thermal

        if actor.state[end] == :L_mode
            # run WPED

             
            α_t = HL_tanh_response(par.tau_t,actor.t_hl, dd.global_time;time_steps=100)
            α_n = HL_tanh_response(par.tau_n,actor.t_hl, dd.global_time;time_steps=100)
            @show alpha_t, alpha_n

            actor.ped_actor = actor.wped_actor
            actor.ped_actor.density_factor = α_n * density_ratio_L_over_H
            run_selected_pedstal_model(actor)
            
            # use cp1d_old as current
            #pedestal = alpha_t .* WPED .+ (1-\alpaha ) .* T_current
        else
            # H mode

            alpha_t = LH_tanh_response(par.tau_t,actor.t_lh, dd.global_time;time_steps=100)
            alpha_n = LH_tanh_response(par.tau_n,actor.t_lh, dd.global_time;time_steps=100)
            @show alpha_t, alpha_n
            
            # t_inf makes alpha = 1
            #t_profile_now = alpha_t .* EPED .+ (1-\alpaha ) .* T_current   
            
            
            #alpha_n = alphan @ t_hl
            # this is going to use density_multiplier

            # pulse schedule will have Hmode density

            # l mode density is then H mode density / multiplier
            # run EPED
            actor.ped_actor = actor.eped_actor
           
            actor.ped_actor.par.density_factor =  (1 - par.density_ratio_L_over_H) * alpha_n + par.density_ratio_L_over_H 
            run_selected_pedstal_model(actor)

#            t_profile_new = alpha_t .* EPED .+ (1-\alpaha ) .* T_current 
            Te_ongoing = (1 .- alpha_t) .* actor.cp1d_transition.electrons.temperature .+  alpha_t .* dd.core_profiles.profiles_1d[].electrons.temperature
            plot(actor.cp1d_transition.electrons.temperature,label="Te old")
            display(plot!(dd.core_profiles.profiles_1d[].electrons.temperature,label="te EPED", ls =:dash))
            @show dd.core_profiles.profiles_1d[].electrons.temperature[1]
            display(plot!(Te_ongoing,label= "Te ongoing", ls =:dash))
            

            plot(actor.cp1d_transition.electrons.density_thermal,label="Te old")
            display(plot!(dd.core_profiles.profiles_1d[].electrons.density_thermal,label="ne EPED", ls =:dash))
            @show dd.core_profiles.profiles_1d[].electrons.density_thermal[1]
#            display(plot!(Te_ongoing,label= "Te ongoing", ls =:dash))

        end
    end

    return actor
end


function LH_tanh_response(τ,t_now, t_LH;time_steps=100)
    t, α = LH_tanh_response(τ, t_LH)
    @show α[1], IMAS.interp1d(t,α).(t_now)
    return maximum([minimum([IMAS.interp1d(t,α).(t_now),1.0]),0.0])
end

function LH_tanh_response(τ,t_LH;time_steps=100)
    t = (0:τ/time_steps:τ) .+ t_LH
    α = tanh.((2pi.*(t .- t_LH .- τ / 4.0)) ./ τ)
    α = ((α .- α[1]) ./ (α[end]-α[1]))
    return t, α
end

function run_selected_pedstal_model(actor)
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