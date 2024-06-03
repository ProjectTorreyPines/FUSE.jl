using CSV
using DataFrames

#= ============ =#
#  ActorCosting  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCosting{T<:Real} <: ParametersActorBuild{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:ARIES, :Sheffield], "-", "Costing model"; default=:ARIES)
    construction_start_year::Entry{Int} = Entry{Int}("year", "Year that plant construction begins"; default=Dates.year(Dates.now()))
    future_inflation_rate::Entry{Real} = Entry{Real}("-", "Predicted average rate of future inflation"; default=0.025)
    plant_lifetime::Entry{Real} = Entry{Real}("year", "Lifetime of the plant"; default=40)
    availability::Entry{Real} = Entry{Real}("-", "Availability fraction of the plant"; default=0.8)
    production_increase::Entry{Real} = Entry{Real}("-", "Factor by which production of ReBCO multiplies"; default=1.0)
    learning_rate::Entry{Real} = Entry{Real}("-", "Learning rate for ReBCO technology production"; default=0.85)
    power_uncertainty::Entry{T} = Entry{T}("-", "Uncertainty to introduce on resulting powers from BOP"; default=0.0)
end

mutable struct ActorCosting{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCosting{P}
    act::ParametersAllActors
    cst_actor::Union{ActorCostingSheffield{D,P},ActorCostingARIES{D,P}}
end

"""
	ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the cost of building, operating, and recommission the fusion power plant.

!!! note 
	Stores data in `dd.costing`
"""
function ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCosting(dd, act.ActorCosting, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCosting(dd::IMAS.dd, par::FUSEparameters__ActorCosting, act::ParametersAllActors; kw...)
    logging_actor_init(ActorCosting)
    par = par(kw...)

    empty!(dd.costing)

    dd.costing.construction_start_year = act.ActorCosting.construction_start_year
    dd.costing.future.inflation_rate = act.ActorCosting.future_inflation_rate
    dd.costing.plant_lifetime = act.ActorCosting.plant_lifetime
    dd.costing.availability = act.ActorCosting.availability
    dd.costing.future.learning.hts.production_increase = act.ActorCosting.production_increase
    dd.costing.future.learning.hts.learning_rate = act.ActorCosting.learning_rate
    dd.costing.power_uncertainty = act.ActorCosting.power_uncertainty

    if par.model == :Sheffield
        cst_actor = ActorCostingSheffield(dd, act.ActorCostingSheffield)
    elseif par.model == :ARIES
        cst_actor = ActorCostingARIES(dd, act.ActorCostingARIES)
    end

    return ActorCosting(dd, par, act, cst_actor)
end

function _step(actor::ActorCosting)
    step(actor.cst_actor)
    return actor
end

function _finalize(actor::ActorCosting)
    store_nominal_values(actor.dd)
    finalize(actor.cst_actor)
    return actor
end
