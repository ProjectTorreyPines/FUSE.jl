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
    model::Switch{Symbol} = Switch{Symbol}([:EPED, :WPED, :auto], "-", "Pedestal model to use"; default=:EPED)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    ne_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_ped_from::Switch{Symbol} = switch_get_from(:zeff_ped)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorPedestal{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPedestal{P}
    ped_actor::Union{Nothing,ActorEPED{D,P},ActorWPED{D,P}}
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
    eped_actor = ActorEPED(dd, act.ActorEPED; ne_ped_from=par.ne_from, par.zeff_ped_from, par.βn_from, par.ip_from, par.rho_nml, par.rho_ped)
    wped_actor = ActorWPED(dd, act.ActorWPED; ne_ped_from=par.ne_from, par.zeff_ped_from, par.rho_ped, par.do_plot)
    return ActorPedestal(dd, par, nothing, eped_actor, wped_actor)
end

"""
    _step(actor::ActorPedestal)

Runs pedestal actor to evaluate pedestal width and height
"""
function _step(actor::ActorPedestal{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    if par.model == :EPED
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

    elseif par.ne_from == :pulse_schedule && par.density_match == :ne_line
        actor.ped_actor.par.ne_ped_from = :core_profiles
        finalize(step(actor.ped_actor))

        summary_ped = dd.summary.local.pedestal
        ne_initial = IMAS.get_from(dd, Val{:ne_ped}, :core_profiles, nothing)
        ne_line_wanted = IMAS.n_e_line(dd.pulse_schedule)
        function cost_ne_ped_from_nel(density_scale, ne_line_wanted)
            @ddtime summary_ped.n_e.value = ne_initial * density_scale

            IMAS.blend_core_edge(mode, cp1d, summary_ped, par.rho_nml, par.rho_ped; what=:densities)

            ne_line = IMAS.geometric_midplane_line_averaged_density(eqt, cp1d)
            return ((ne_line - ne_line_wanted) / ne_line_wanted)^2
        end
        res = Optim.optimize(x -> cost_ne_ped_from_nel(x, ne_line_wanted), 0.01, 100, Optim.GoldenSection(); rel_tol=1E-3)
        cost_ne_ped_from_nel(res.minimizer, ne_line_wanted)

        finalize(step(actor.ped_actor))
        actor.ped_actor.par.ne_ped_from = :pulse_schedule
    end

    return actor
end
