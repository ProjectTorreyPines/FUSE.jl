using CSV
using DataFrames

#= ============ =#
#  ActorCosting  #
#= ============ =#
@actor_parameters_struct ActorCosting{T} begin
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
    par::OverrideParameters{P,FUSEparameters__ActorCosting{P}}
    act::ParametersAllActors{P}
    cst_actor::Union{ActorCostingSheffield{Measurement{Float64},P},ActorCostingARIES{Measurement{Float64},P}}
end

"""
    ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the comprehensive cost of building, operating, and decommissioning a fusion power plant.

This actor provides a unified interface to different costing methodologies and handles:
- Direct capital costs for all major plant systems and components  
- Operating and maintenance costs over the plant lifetime
- Decommissioning costs at end-of-life
- Economic parameters like inflation, learning rates, and financial factors

Available costing models:
- `:ARIES`: Based on ARIES tokamak reactor studies with detailed component breakdown
- `:Sheffield`: Based on Sheffield & Milora generic magnetic fusion reactor methodology

The actor creates a separate data structure with uncertainty propagation (Measurement types)
to handle cost uncertainties, then transfers final results back to the main data structure.
Key economic outputs include levelized cost of electricity (LCOE) and lifetime plant costs.

!!! note

    Results are stored in `dd.costing` with detailed breakdown by system and subsystem
"""
function ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCosting(dd, act.ActorCosting, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCosting(dd::IMAS.dd, par::FUSEparameters__ActorCosting, act::ParametersAllActors; kw...)
    logging_actor_init(ActorCosting)
    par = OverrideParameters(par; kw...)

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
