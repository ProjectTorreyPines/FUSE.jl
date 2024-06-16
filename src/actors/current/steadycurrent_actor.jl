#= ======================= =#
#  ActorSteadyStateCurrent  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorSteadyStateCurrent{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    allow_floating_plasma_current::Entry{Bool} = Entry{Bool}("-", "Zero loop voltage if non-inductive fraction exceeds 100% of the target Ip")
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
end

mutable struct ActorSteadyStateCurrent{D,P} <: SingleAbstractActor{D,P}
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
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cpg = dd.core_profiles.global_quantities
    cp1d = dd.core_profiles.profiles_1d[]

    ip_target = IMAS.get_from(dd, Val{:ip}, par.ip_from)

    # update j_ohmic
    IMAS.j_ohmic_steady_state!(eqt, dd.core_profiles.profiles_1d[], ip_target)

    # allow floating plasma current
    if abs(ip_target) < abs(@ddtime(cpg.current_non_inductive)) && par.allow_floating_plasma_current
        cp1d.j_ohmic = zeros(length(cp1d.grid.rho_tor_norm))
    end

    return actor
end
