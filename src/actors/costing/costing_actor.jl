using CSV
using DataFrames
import Memoize

#= ============ =#
#  ActorCosting  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCosting{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:FUSE, :ARIES, :Sheffield], "-", "Costing model"; default=:ARIES)
end

mutable struct ActorCosting <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorCosting
    act::ParametersAllActors
    cst_actor::Union{ActorSheffieldCosting,ActorARIESCosting}
end

"""
    ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the cost of building, operating, and recommission the fusion power plant.

!!! note 
    Stores data in `dd.costing`
"""
function ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCosting(kw...)
    actor = ActorCosting(dd, par, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCosting(dd::IMAS.dd, par::FUSEparameters__ActorCosting, act::ParametersAllActors; kw...)
    logging_actor_init(ActorCosting)
    par = par(kw...)
    if par.model == :Sheffield
        cst_actor = ActorSheffieldCosting(dd, act.ActorSheffieldCosting)
    elseif par.model == :ARIES
        cst_actor = ActorARIESCosting(dd, act.ActorARIESCosting)
    else 
        error("ActorCosting: model = '$(par.model)' can only be ':Sheffield' or ':ARIES'")
    end
    return ActorCosting(dd, par, act, cst_actor)
end

function _step(actor::ActorCosting)
    step(actor.cst_actor)
    return actor
end

function _finalize(actor::ActorCosting)
    finalize(actor.cst_actor)
    return actor
end
