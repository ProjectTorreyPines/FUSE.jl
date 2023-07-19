#= ================== =#
#  ActorWholeFacility  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorWholeFacility{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    update_plasma::Entry{Bool} = Entry{Bool}("-", "Run plasma related actors"; default=true)
end

mutable struct ActorWholeFacility{D,P} <: FacilityAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorWholeFacility{P}
    act::ParametersAllActors
    EquilibriumTransport::Union{Nothing,ActorEquilibriumTransport{D,P}}
    StabilityLimits::Union{Nothing,ActorStabilityLimits{D,P}}
    HFSsizing::Union{Nothing,ActorHFSsizing{D,P}}
    LFSsizing::Union{Nothing,ActorLFSsizing{D,P}}
    CXbuild::Union{Nothing,ActorCXbuild{D,P}}
    PFcoilsOpt::Union{Nothing,ActorPFcoilsOpt{D,P}}
    PassiveStructures::Union{Nothing,ActorPassiveStructures{D,P}}
    Neutronics::Union{Nothing,ActorNeutronics{D,P}}
    Blanket::Union{Nothing,ActorBlanket{D,P}}
    Divertors::Union{Nothing,ActorDivertors{D,P}}
    BalanceOfPlant::Union{Nothing,ActorBalanceOfPlant{D,P}}
    Costing::Union{Nothing,ActorCosting{D,P}}
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
    actor = ActorWholeFacility(dd, act.ActorWholeFacility, act; kw...)
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
        actor.StabilityLimits = ActorStabilityLimits(dd, act)
    end

    actor.HFSsizing = ActorHFSsizing(dd, act)
    if abs(actor.HFSsizing.R0_scale - 1.0) > 1E-6
        if !par.update_plasma
            error("HFSsizing has changed aspect ratio by $((actor.HFSsizing.R0_scale-1.0)*100)%. You must allow for `act.ActorWholeFacility.update_plasma=true`.")
        end
        @warn "Aspect ratio changed by $((actor.HFSsizing.R0_scale-1.0)*100)% --> re-running plasma actors"
        scale_aspect_ratio!(dd, actor.HFSsizing.R0_scale)
        actor.EquilibriumTransport = ActorEquilibriumTransport(dd, act)
        actor.StabilityLimits = ActorStabilityLimits(dd, act)
    end

    actor.LFSsizing = ActorLFSsizing(dd, act)

    actor.CXbuild = ActorCXbuild(dd, act)

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

"""
    scale_aspect_ratio!(ps::IMAS.pulse_schedule, R0_scale::Float64)

Update the radial coordinates of all the leaves of pulse_schedule and wall IDSs
as well as equilibrium.vacuum_toroidal_field.r0
"""
function scale_aspect_ratio!(dd::IMAS.dd, R0_scale::Float64)
    eq = dd.equilibrium
    ΔR0 = eq.vacuum_toroidal_field.r0 * (R0_scale - 1.0)
    eq.vacuum_toroidal_field.r0 += ΔR0

    dd.pulse_schedule.tf.b_field_tor_vacuum_r.reference.data *= R0_scale

    for leaf in (collect(IMAS.leaves(dd.pulse_schedule)); collect(IMAS.leaves(dd.wall)))
        uloc = IMAS.ulocation(leaf.ids, leaf.field)
        if occursin(".r.", uloc) && IMAS.info(uloc)["units"] ∈ ("mixed", "m")
            old_value = getproperty(leaf.ids, leaf.field)
            setproperty!(leaf.ids, leaf.field, old_value .+ ΔR0)
        end
    end
end