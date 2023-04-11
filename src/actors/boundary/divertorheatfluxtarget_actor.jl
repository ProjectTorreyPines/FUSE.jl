#= ========================= =#
#  ActorDivertorHeatFluxTarget #
#= ========================= =#

import BoundaryPlasmaModels
import BoundaryPlasmaModels: FUSEparameters__ActorDivertorHeatFluxTarget

#==========#

mutable struct ActorDivertorHeatFluxTarget <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorDivertorHeatFluxTarget
    model::Union{Nothing,BoundaryPlasmaModels.DivertorHeatFluxModel}
end

"""
ActorDivertorHeatFluxTarget(dd::IMAS.dd, act::ParametersAllActors; kw...)

doc
"""
function ActorDivertorHeatFluxTarget(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorDivertorHeatFluxTarget(kw...)
    actor = ActorDivertorHeatFluxTarget(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorDivertorHeatFluxTarget(dd::IMAS.dd, par::FUSEparameters__ActorDivertorHeatFluxTarget; kw...)
    logging_actor_init(ActorDivertorHeatFluxTarget)
    par = par(kw...)
    return ActorDivertorHeatFluxTarget(dd, par, nothing)
end

function _step(actor::ActorDivertorHeatFluxTarget)
    actor.model = BoundaryPlasmaModels.DivertorHeatFluxModel(actor.dd,actor.par)
    return actor
end

function _finalize(actor::ActorDivertorHeatFluxTarget; kw...)
    BoundaryPlasmaModels.export2dd(actor.dd, actor.model; kw...)
    actor
end