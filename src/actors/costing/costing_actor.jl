using CSV
using DataFrames
using Memoize: Memoize

#= ============ =#
#  ActorCosting  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCosting{T} <: ParametersActor where {T <: Real}
	_parent::WeakRef = WeakRef(nothing)
	_name::Symbol = :not_set
	model::Switch{Symbol} = Switch(Symbol, [:FUSE, :ARIES, :Sheffield], "-", "Costing model"; default = :ARIES)
	construction_start_year::Entry{Int} = Entry(Int, "year", "Year that plant construction begins"; default = Dates.year(Dates.now()))
	future_inflation_rate::Entry{T} = Entry(T, "-", "Predicted average rate of future inflation"; default = 0.025)
	plant_lifetime::Entry{Int} = Entry(Int, "year", "Lifetime of the plant"; default = 40)
	availability::Entry{T} = Entry(T, "-", "Availability fraction of the plant"; default = 0.8)
	production_increase::Entry{T} = Entry(T, "-", "Factor by which production of ReBCO multiplies"; default = 10)
	learning_rate::Entry{T} = Entry(T, "-", "Learning rate for ReBCO technology production"; default = 0.85)
end

mutable struct ActorCosting <: FacilityAbstractActor
	dd::IMAS.dd
	par::FUSEparameters__ActorCosting
	act::ParametersAllActors
	cst_actor::Union{ActorSheffieldCosting, ActorARIESCosting}
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
