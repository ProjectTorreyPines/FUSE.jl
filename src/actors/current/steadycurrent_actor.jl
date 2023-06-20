#= ======================= =#
#  ActorSteadyStateCurrent  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorSteadyStateCurrent{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    allow_floating_plasma_current::Entry{Bool} = Entry{Bool}("-", "allows the plasma current to increase or decrease based on the non-inductive current"; default=false)
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

Evolves the current to steady state using the conductivity from `dd.core_profiles` and current profile form `dd.equilibrium`.

Also sets the ohmic, bootstrap and non-inductive current profiles in `dd.core_profiles`

!!! note 
    Stores data in `dd.core_profiles`, `dd.equilbrium`
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

    # update j_ohmic (this also restores j_tor, j_total as expressions)
    IMAS.j_ohmic_steady_state!(eqt, dd.core_profiles.profiles_1d[])
    # update core_sources related to current

    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)
    return actor
end
