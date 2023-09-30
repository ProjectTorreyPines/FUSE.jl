#= ======================= =#
#  ActorSteadyStateCurrent  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorSteadyStateCurrent{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    allow_floating_plasma_current::Entry{Bool} = Entry{Bool}("-", "allows the plasma current to increase or decrease based on the non-inductive current"; default=false)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
end

mutable struct ActorSteadyStateCurrent{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorSteadyStateCurrent{P}
end

"""
    ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

* evolves the ohmic current to steady state using the conductivity from `dd.core_profiles` and total current form `dd.equilibrium`.
* sets the ohmic, bootstrap, and non-inductive current profiles in `dd.core_profiles`
* updates bootstrap and ohmic in `dd.core_sources`

!!! note 
    Stores data in `dd.core_sources` and `dd.core_profiles`
"""
function ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorSteadyStateCurrent(dd, par::FUSEparameters__ActorSteadyStateCurrent; kw...)
    logging_actor_init(ActorSteadyStateCurrent)
    par = par(kw...)
    return ActorSteadyStateCurrent(dd, par)
end

function _step(actor::ActorSteadyStateCurrent)
    dd = actor.dd
    par = actor.par
    eqt = dd.equilibrium.time_slice[]
    cpg = dd.core_profiles.global_quantities
    cp1d = dd.core_profiles.profiles_1d[]

    ip_target = IMAS.get_from(dd, Val{:ip}, par.ip_from)
    if abs(ip_target) < abs(@ddtime(cpg.current_non_inductive))
        if par.allow_floating_plasma_current
            println("set j_ohmic to zero and allow ip to be floating")
            cp1d.j_ohmic = zeros(length(cp1d.grid.rho_tor_norm))
        else
            @warn "j_ohmic will change sign, as the non-inductive current $(round(@ddtime(cpg.current_non_inductive),digits=3)) has larger magnitude than ip_target $(round(ip_target,digits=3))"
            IMAS.j_ohmic_steady_state!(eqt, dd.core_profiles.profiles_1d[], ip_target)
        end
    else
        # update j_ohmic (this also restores j_tor, j_total as expressions)
        IMAS.j_ohmic_steady_state!(eqt, dd.core_profiles.profiles_1d[], ip_target)
    end

    # update core_sources related to current
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)

    # add vloop info to pulse_schedule
    vloop = IMAS.get_from(dd, Val{:vloop}, :core_profiles)
    @ddtime(dd.pulse_schedule.flux_control.loop_voltage.reference.data = vloop)

    return actor
end
