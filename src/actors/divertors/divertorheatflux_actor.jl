#= ===================== =#
#  ActorDivertorHeatFlux  #
#= ===================== =#

import BoundaryPlasmaModels

Base.@kwdef mutable struct FUSEparameters__ActorDivertorHeatFlux{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:lengyel], "-", "Divertor heat flux model"; default=:lengyel)
end

mutable struct ActorDivertorHeatFlux <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorDivertorHeatFlux
    boundary_plasma_model::Union{BoundaryPlasmaModels.DivertorHeatFluxModel,Nothing}
end

"""
    ActorDivertorHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates divertor heat flux

!!! note 
    Manipulates data in ???
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
    dd = actor.dd
    par = actor.par
    actor.boundary_plasma_model = BoundaryPlasmaModels.DivertorHeatFluxModel(par.model) #create instance of divertor heat flux model 
    BoundaryPlasmaModels.setup_model(actor.boundary_plasma_model, dd) # load model parameters from dd into dhf
    actor.boundary_plasma_model() #execute model
    return actor
end

function _finalize(actor::ActorDivertorHeatFlux)
    dd = actor.dd
    BoundaryPlasmaModels.export2dd(dd, actor.boundary_plasma_model) # export to dd structures
    return actor
end