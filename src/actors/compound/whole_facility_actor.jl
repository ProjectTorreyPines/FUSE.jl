#= ================== =#
#  ActorWholeFacility  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorWholeFacility{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    update_plasma::Entry{Bool} = Entry(Bool, "-", "Run plasma related actors"; default=true)
end

mutable struct ActorWholeFacility <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorWholeFacility
    act::ParametersAllActors
    EquilibriumTransport::Union{Nothing,ActorEquilibriumTransport}
    StabilityLimits::Union{Nothing,ActorStabilityLimits}
    HFSsizing::Union{Nothing,ActorHFSsizing}
    LFSsizing::Union{Nothing,ActorLFSsizing}
    CXbuild::Union{Nothing,ActorCXbuild}
    PFcoilsOpt::Union{Nothing,ActorPFcoilsOpt}
    PassiveStructures::Union{Nothing,ActorPassiveStructures}
    Neutronics::Union{Nothing,ActorNeutronics}
    Blanket::Union{Nothing,ActorBlanket}
    Divertors::Union{Nothing,ActorDivertors}
    BalanceOfPlant::Union{Nothing,ActorBalanceOfPlant}
    Costing::Union{Nothing,ActorCosting}
end

"""
    ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs all the physics, engineering and costing actors needed to model the whole plant:
* ActorEquilibriumTransport
    * ActorSteadyStateCurrent
    * ActorHCD
    * ActorCoreTransport
    * ActorEquilibrium
* ActorStabilityLimits
* ActorHFSsizing
* ActorLFSsizing
* ActorCXbuild
* ActorPFcoilsOpt
* ActorPassiveStructures
* ActorNeutronics
* ActorBlanket
* ActorDivertors
* ActorBalanceOfPlant
* ActorCosting

!!! note 
    Stores data in `dd`
"""
function ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorWholeFacility
    actor = ActorWholeFacility(dd, par, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorWholeFacility(dd::IMAS.dd, par::FUSEparameters__ActorWholeFacility, act::ParametersAllActors; kw...)
    logging_actor_init(ActorWholeFacility)
    par = par(kw...)

    ActorWholeFacility(dd, par, act,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing)
end

function _step(actor::ActorWholeFacility)
    dd = actor.dd
    par = actor.par
    act = actor.act
    if par.update_plasma
        actor.EquilibriumTransport = ActorEquilibriumTransport(dd, act)
        actor.StabilityLimits == ActorStabilityLimits(dd, act)
    end
    actor.HFSsizing = ActorHFSsizing(dd, act)
    actor.LFSsizing = ActorLFSsizing(dd, act)
    actor.CXbuild = ActorCXbuild(dd, act)
    if par.update_plasma
        # depending on act.ActorHFSsizing.aspect_ratio_tolerance setting the wall may shift slightly with respect to the equilibrium
        ActorEquilibrium(dd, act)
    end
    actor.PFcoilsOpt = ActorPFcoilsOpt(dd, act)
    if act.ActorPFcoilsOpt.update_equilibrium && act.ActorCXbuild.rebuild_wall
        actor.CXbuild = ActorCXbuild(dd, act)
        actor.PFcoilsOpt = ActorPFcoilsOpt(dd, act; update_equilibrium=false)
    end
    actor.PassiveStructures = ActorPassiveStructures(dd, act)
    actor.Neutronics = ActorNeutronics(dd, act)
    actor.Blanket = ActorBlanket(dd, act)
    actor.Divertors = ActorDivertors(dd, act)
    actor.BalanceOfPlant = ActorBalanceOfPlant(dd, act)
    actor.Costing = ActorCosting(dd, act)
    return actor
end
