#= ======================= =#
#  ActorSteadyStateCurrent  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorSteadyStateCurrent{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorSteadyStateCurrent{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorSteadyStateCurrent{P}
    function ActorSteadyStateCurrent(dd::IMAS.dd{D}, par::FUSEparameters__ActorSteadyStateCurrent{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSteadyStateCurrent)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
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

function _step(actor::ActorSteadyStateCurrent)
    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    # update j_ohmic (this also restores j_tor, j_total as expressions)
    Ip = IMAS.get_time_array(dd.pulse_schedule.flux_control.i_plasma.reference, :data, dd.global_time, :linear)
    IMAS.j_ohmic_steady_state!(eqt, cp1d)#, Ip)
    # update core_sources related to current
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)
    return actor
end
