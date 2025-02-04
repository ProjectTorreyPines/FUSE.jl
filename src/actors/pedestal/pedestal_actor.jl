import EPEDNN

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
    model::Switch{Symbol} = Switch{Symbol}([:EPED, :WPED, :auto, :none], "-", "Pedestal model to use"; default=:EPED)
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
    return ActorPedestal(dd, par, act, none_actor, none_actor, eped_actor, wped_actor)
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
        mode = :H_mode
    elseif par.model == :WPED
        actor.ped_actor = actor.wped_actor
        mode = :L_mode
    elseif par.model == :auto
        if IMAS.satisfies_h_mode_conditions(dd)
            actor.ped_actor = actor.eped_actor
            mode = :H_mode
        else
            actor.ped_actor = actor.wped_actor
            mode = :L_mode
            if eqt.boundary.triangularity < 0.0
                actor.wped_actor.par.ped_to_core_fraction = 0.3
            else
                actor.wped_actor.par.ped_to_core_fraction = 0.05
            end
        end
    end

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

        cp1d_copy = deepcopy(cp1d)
        factor = 1.0
        for k in 1:10
            # scale thermal densities to match desired line average (and temperatures accordingly, in case they matter)
            ne_line = IMAS.geometric_midplane_line_averaged_density(eqt, cp1d)
            ne_line_wanted = IMAS.ne_line(dd.pulse_schedule)
            factor = factor * ne_line_wanted / ne_line
            cp1d.electrons.density_thermal = cp1d_copy.electrons.density_thermal * factor
            for i in eachindex(cp1d.ion)
                cp1d.ion[i].density_thermal = cp1d_copy.ion[i].density_thermal * factor
            end
            cp1d.electrons.temperature = cp1d_copy.electrons.temperature / factor
            for i in eachindex(cp1d.ion)
                cp1d.ion[i].temperature = cp1d_copy.ion[i].temperature / factor
            end

            # run the pedestal model
            if k == 1
                finalize(step(actor.ped_actor))
            else
                _finalize(_step(actor.ped_actor))
            end

            ne_line = IMAS.geometric_midplane_line_averaged_density(eqt, cp1d)
            ne_line_wanted = IMAS.ne_line(dd.pulse_schedule)
            if abs(ne_line - ne_line_wanted) / ne_line_wanted < 1E-3
                break
            end
        end

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