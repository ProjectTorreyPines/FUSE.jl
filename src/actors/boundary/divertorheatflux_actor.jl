#= ========================= =#
#  ActorDivertorHeatFlux      #
#= ========================= =#

import BoundaryPlasmaModels
import BoundaryPlasmaModels: FUSEparameters__ActorDivertorHeatFlux


#==========#

mutable struct ActorDivertorHeatFlux <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorDivertorHeatFlux
    model::Union{BoundaryPlasmaModels.DivertorHeatFluxModel,Nothing}
end

"""
ActorDivertorHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)

doc
"""
function ActorDivertorHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorDivertorHeatFlux(kw...)
    actor = ActorDivertorHeatFlux(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorDivertorHeatFlux(dd::IMAS.dd, par::FUSEparameters__ActorDivertorHeatFlux; kw...)
    logging_actor_init(ActorDivertorHeatFlux)
    par = par(kw...)
    return ActorDivertorHeatFlux(dd, par, nothing)
end

function _step(actor::ActorDivertorHeatFlux)
    actor.model = BoundaryPlasmaModels.DivertorHeatFluxModel(actor.par) #create instance of divertor heat flux model 
    BoundaryPlasmaModels.setup_model(actor.model,actor.dd) # load model parameters from dd into dhf
    actor.model() #execute model
    return actor
end

function _finalize(actor::ActorDivertorHeatFlux; kw...)
    BoundaryPlasmaModels.export2dd(actor.dd, actor.model; kw...) # export to dd structures
    actor
end