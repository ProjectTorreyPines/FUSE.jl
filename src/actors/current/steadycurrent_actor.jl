#= ======================= =#
#  ActorSteadyStateCurrent  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorSteadyStateCurrent{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorSteadyStateCurrent <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorSteadyStateCurrent
    function ActorSteadyStateCurrent(dd::IMAS.dd, par::FUSEparameters__ActorSteadyStateCurrent; kw...)
        logging_actor_init(ActorSteadyStateCurrent)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the current to steady state using the conductivity from `dd.core_profiles` and current profile form `dd.equilibrium`.

Also sets the ohmic, bootstrap and non-inductive current profiles in `dd.core_profiles`

!!! note 
    Stores data in `dd.core_profiles, dd.equilbrium`
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
    IMAS.j_ohmic_steady_state!(eqt, cp1d)
    # update core_sources related to current
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)
    return actor
end
