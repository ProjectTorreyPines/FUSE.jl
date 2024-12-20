#= ================== =#
#  ActorWholeFacility  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorWholeFacility{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    update_plasma::Entry{Bool} = Entry{Bool}("-", "Run plasma related actors"; default=true)
    update_build::Entry{Bool} = Entry{Bool}("-", "Optimize tokamak build"; default=true)
end

mutable struct ActorWholeFacility{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorWholeFacility{P}
    act::ParametersAllActors
    StationaryPlasma::Union{Nothing,ActorStationaryPlasma{D,P}}
    StabilityLimits::Union{Nothing,ActorStabilityLimits{D,P}}
    FluxSwing::Union{Nothing,ActorFluxSwing{D,P}}
    Stresses::Union{Nothing,ActorStresses{D,P}}
    HFSsizing::Union{Nothing,ActorHFSsizing{D,P}}
    LFSsizing::Union{Nothing,ActorLFSsizing{D,P}}
    CXbuild::Union{Nothing,ActorCXbuild{D,P}}
    PFdesign::Union{Nothing,ActorPFdesign{D,P}}
    PFactive::Union{Nothing,ActorPFactive{D,P}}
    PassiveStructures::Union{Nothing,ActorPassiveStructures{D,P}}
    VerticalStability::Union{Nothing,ActorVerticalStability{D,P}}
    Neutronics::Union{Nothing,ActorNeutronics{D,P}}
    Blanket::Union{Nothing,ActorBlanket{D,P}}
    Divertors::Union{Nothing,ActorDivertors{D,P}}
    BalanceOfPlant::Union{Nothing,ActorBalanceOfPlant{D,P}}
    Costing::Union{Nothing,ActorCosting{D,P}}
end

"""
    ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs all the physics, engineering and costing actors needed to model the whole plant:

  - ActorStationaryPlasma
  - ActorStabilityLimits
  - ActorHFSsizing
  - ActorLFSsizing
  - ActorCXbuild
  - ActorFluxSwing
  - ActorStresses
  - ActorPFdesign
  - ActorPFactive
  - ActorPassiveStructures
  - ActorVerticalStability
  - ActorNeutronics
  - ActorBlanket
  - ActorDivertors
  - ActorBalanceOfPlant
  - ActorCosting

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

    return ActorWholeFacility(dd, par, act,
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
        nothing,
        nothing,
        nothing,
        nothing,
        nothing)
end

function _step(actor::ActorWholeFacility)
    dd = actor.dd
    par = actor.par
    act = deepcopy(actor.act)

    if par.update_build && act.ActorCXbuild.rebuild_wall
        # when we rebuild the wall, we let the strike points do whatever they want, since the wall will adapt instead
        act.ActorPFactive.strike_points_weight = 0.0
    end

    if !isempty(dd.build.layer) && par.update_build
        # we start by optimizing coil location, so that when we go solve the equilibrium we can hold it in place
        actor.PFdesign = ActorPFdesign(dd, act)
    end

    if par.update_plasma
        # remove the wall from the calculation to avoid switching between limited/diverted plasma
        wall_bkp = deepcopy(dd.wall)
        empty!(dd.wall)
        if !isempty(dd.pf_active.coil)
            IMAS.first_wall!(dd.wall, IMAS.first_wall(dd.pf_active)...)
        end
        actor.StationaryPlasma = ActorStationaryPlasma(dd, act)
        dd.wall = wall_bkp
    end

    if isempty(dd.build.layer)
        @warn "ActorWholeFacility: skipping engineering/costing actors since build is missing"

    else
        if par.update_build
            actor.HFSsizing = ActorHFSsizing(dd, act)
            actor.LFSsizing = ActorLFSsizing(dd, act)
            empty!(dd.pf_active) # remove old coils since those affect the shape of the layers
            actor.CXbuild = ActorCXbuild(dd, act)

            actor.PFdesign = ActorPFdesign(dd, act)
            if act.ActorCXbuild.rebuild_wall
                ActorEquilibrium(dd, act; ip_from=:equilibrium, j_p_from=:equilibrium)
                actor.CXbuild = ActorCXbuild(dd, act)
            end

        else
            actor.FluxSwing = ActorFluxSwing(dd, act)
            actor.Stresses = ActorStresses(dd, act)
            actor.PFactive = ActorPFactive(dd, act)
        end

        actor.Neutronics = ActorNeutronics(dd, act)

        actor.Blanket = ActorBlanket(dd, act)
        actor.CXbuild = ActorCXbuild(dd, act)

        actor.PassiveStructures = ActorPassiveStructures(dd, act)

        actor.Divertors = ActorDivertors(dd, act)

        actor.StabilityLimits = ActorStabilityLimits(dd, act)
        actor.VerticalStability = ActorVerticalStability(dd, act)

        actor.BalanceOfPlant = ActorBalanceOfPlant(dd, act)

        actor.Costing = ActorCosting(dd, act)
    end

    return actor
end
