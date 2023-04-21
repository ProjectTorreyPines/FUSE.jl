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
    construction_start_year::Entry{T} = Entry(T, "-", "Year that plant construction begins"; default=2030)
    construction_lead_time::Entry{T} = Entry(T, "years", "Duration of construction"; default=8)
    inflate_to_start_year::Entry{Bool} = Entry(Bool, "-", "Return costs in dollars inflated to year that construction begins"; default = false)
    future_inflation_rate::Entry{T} = Entry(T, "-", "Predicted average rate of future inflation"; default = 0.025)
    land_space::Entry{T} = Entry(T, "acres", "Plant site space required in acres"; default=1000.0)
    building_volume::Entry{T} = Entry(T, "m^3", "Volume of the tokmak building"; default=140.0e3)
    interest_rate::Entry{T} = Entry(T, "-", "Annual interest rate fraction of direct capital cost"; default=0.05)
    fixed_charge_rate::Entry{T} = Entry(T, "-", "Constant dollar fixed charge rate"; default=0.078)
    indirect_cost_rate::Entry{T} = Entry(T, "-", "Indirect cost associated with construction, equipment, services, engineering construction management and owners cost"; default=0.4) 
    lifetime::Entry{Int} = Entry(Int, "years", "lifetime of the plant"; default=40)
    availability::Entry{T} = Entry(T, "-", "availability fraction of the plant"; default=0.803)
    escalation_fraction::Entry{T} = Entry(T, "-", "yearly escalation fraction based on risk assessment"; default=0.05)
    blanket_lifetime::Entry{T} = Entry(T, "years", "lifetime of the blanket"; default=6.8)
    initial_cost_blanket::Entry{T} = Entry(T, "millions of dollars", "cost of initial blanket"; default=200)
    initial_cost_divertor::Entry{T} = Entry(T, "millions of dollars", "cost of initial divertor"; default=8)
    divertor_fluence_lifetime::Entry{T} = Entry(T, "MW/yr*m^2", "divertor fluence lifetime"; default=10)
    blanket_fluence_lifetime::Entry{T} = Entry(T, "MW/yr*m^2", "blanket fluence lifetime"; default=15)
end

mutable struct ActorCosting <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorCosting
    function ActorCosting(dd::IMAS.dd, par::FUSEparameters__ActorCosting; kw...)
        logging_actor_init(ActorCosting)
        par = par(kw...)
        return new(dd, par)
    end
end


#= ==================== =#
#  Inflation Adjustment  #
#= ==================== =#

mutable struct DollarAdjust
    future_inflation_rate::Real 
    construction_start_year::Int
    inflate_to_start_year::Bool 
    year_assessed::Union{Missing,Int}
    year::Union{Missing,Int}
end

function DollarAdjust(par::FUSEparameters__ActorCosting)
    return DollarAdjust(par.future_inflation_rate, par.construction_start_year, par.inflate_to_start_year, missing, missing)
end 

Memoize.@memoize function load_inflation_rate()
    csv_loc = abspath(joinpath(@__DIR__, "CPI.csv"))
    CPI = DataFrame(CSV.File(csv_loc))
    return CPI
end

"""
    future_dollars(dollars::Real, da::DollarAdjust)

    Adjusts costs assessed in a previous year to current or future dollar amount 

    NOTE: Inflation of past costs according to U.S. Bureau of Labor Statistics CPI data 
    Source: https://data.bls.gov/cgi-bin/surveymost?cu (Dataset CUUR0000SA0)
"""

function future_dollars(dollars::Real, da::DollarAdjust)
    CPI = load_inflation_rate()
    if 1913 <= da.year_assessed <= 2022 
        index = CPI.Year .== da.year_assessed 
        CPI_past_year = CPI[index, "Year Avg"][1]
        index_2023 = CPI.Year .== 2023 
        val_today = CPI[index_2023, "Year Avg"][1] ./ CPI_past_year .* dollars
    elseif da.year_assessed == 2023
        val_today = dollars 
    else 
        @warn "Inflation data not available for $(da.year_assessed)"
        val_today = dollars 
    end

    if da.inflate_to_start_year
        if da.construction_start_year < 2023
            @warn "inflate_to_start_year = true requires that construction_start_year is in the future. Costs will be returned in 2023 dollars."
            year = 2023
        else
            year = da.construction_start_year 
        end
    else 
        year = 2023
    end

    n_years = year - 2023
    future_val = val_today * ((1 + da.future_inflation_rate)^n_years)

    da.year_assessed = missing #wipe out year_assessed each time to avoid propogating the wrong year 

    return future_val
end

"""
    ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the cost of building, operating, and recommission the fusion power plant.

!!! note 
    Stores data in `dd.costing`
"""
function ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCosting
    actor = ActorCosting(dd, par; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorCosting)
    par = actor.par
    dd = actor.dd
    cst = dd.costing

    empty!(cst)

    if par.model == :ARIES
        costing_ARIES(dd, par)
    elseif par.model == :Sheffield
        costing_Sheffield(dd, par)
    elseif par.model == :FUSE
        costing_FUSE(dd, par)
    end

    return actor
end

function _finalize(actor::ActorCosting)
    # sort system/subsystem by their costs
    sort!(actor.dd.costing.cost_direct_capital.system, by=x -> x.cost, rev=true)
    for sys in actor.dd.costing.cost_direct_capital.system
        sort!(sys.subsystem, by=x -> x.cost, rev=true)
    end
    return actor
end
