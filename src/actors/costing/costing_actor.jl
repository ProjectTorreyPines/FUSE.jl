using CSV
using DataFrames

#= ============ =#
#  ActorCosting  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCosting{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:ARIES, :Sheffield], "-", "Costing model"; default=:ARIES)
    construction_start_year::Entry{Measurement{Float64}} = Entry{Measurement{Float64}}("year", "Year that plant construction begins"; default=Dates.year(Dates.now()) + 10.0 ± 5.0)
    future_inflation_rate::Entry{Measurement{Float64}} = Entry{Measurement{Float64}}("-", "Predicted average rate of future inflation"; default=0.025 ± 0.0125)
    plant_lifetime::Entry{Measurement{Float64}} = Entry{Measurement{Float64}}("year", "Lifetime of the plant"; default=35.0 ± 5.0)
    availability::Entry{Measurement{Float64}} = Entry{Measurement{Float64}}("-", "Availability fraction of the plant"; default=0.8 ± 0.1)
    production_increase::Entry{Measurement{Float64}} = Entry{Measurement{Float64}}(
        "-",
        "Factor by which production of ReBCO multiplies per year";
        default=1.0 ± 0.5,
        check=x -> @assert x >= 0 "production_increase must be >= 0.0"
    )
    learning_rate::Entry{Measurement{Float64}} = Entry{Measurement{Float64}}(
        "-",
        "Learning rate for ReBCO technology production";
        default=0.85 ± 0.1,
        check=x -> @assert 1.0 >= x >= 0.0 "learning_rate must be between 0.0 and 1.0"
    )
end

mutable struct ActorCosting{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCosting{P}
    act::ParametersAllActors{P}
    cst_actor::Union{ActorCostingSheffield{Measurement{Float64},P},ActorCostingARIES{Measurement{Float64},P}}
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

    ddR = IMAS.get_timeslice(Measurement{Float64}, dd)
    copy_workflow!(ddR, dd) # to let FUSE continue tracking the execution workflow on ddR

    ddR.costing.construction_start_year = act.ActorCosting.construction_start_year
    ddR.costing.future.inflation_rate = act.ActorCosting.future_inflation_rate
    ddR.costing.plant_lifetime = act.ActorCosting.plant_lifetime
    ddR.costing.availability = act.ActorCosting.availability
    ddR.costing.future.learning.hts.production_increase = act.ActorCosting.production_increase
    ddR.costing.future.learning.hts.learning_rate = act.ActorCosting.learning_rate

    if par.model == :Sheffield
        cst_actor = ActorCostingSheffield(ddR, act.ActorCostingSheffield)
    elseif par.model == :ARIES
        cst_actor = ActorCostingARIES(ddR, act.ActorCostingARIES)
    end

    return ActorCosting(dd, par, act, cst_actor)
end

function _step(actor::ActorCosting)
    step(actor.cst_actor)
    return actor
end

function _finalize(actor::ActorCosting)
    finalize(actor.cst_actor)

    cst = actor.cst_actor.dd.costing
    fill!(actor.dd.costing, cst)

    return actor
end
