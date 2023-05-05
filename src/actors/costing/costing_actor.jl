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
    construction_start_year::Entry{Int} = Entry(Int, "year", "Year that plant construction begins"; default=Dates.year(Dates.now()))
    future_inflation_rate::Entry{T} = Entry(T, "-", "Predicted average rate of future inflation"; default=0.025)
    plant_lifetime::Entry{Int} = Entry(Int, "year", "Lifetime of the plant"; default=40)
    availability::Entry{T} = Entry(T, "-", "Availability fraction of the plant"; default=0.8)

    # ARIES
    land_space::Entry{T} = Entry(T, "acres", "Plant site space required"; default=1000.0)
    building_volume::Entry{T} = Entry(T, "m^3", "Volume of the tokmak building"; default=140.0e3)
    interest_rate::Entry{T} = Entry(T, "-", "Annual interest rate fraction of direct capital cost"; default=0.05)
    indirect_cost_rate::Entry{T} = Entry(T, "-", "Indirect cost associated with construction, equipment, services, engineering construction management and owners cost"; default=0.4)
    escalation_fraction::Entry{T} = Entry(T, "-", "Yearly escalation fraction based on risk assessment"; default=0.05)
    blanket_lifetime::Entry{T} = Entry(T, "year", "Lifetime of the blanket"; default=6.8)

    # Sheffield
    construction_lead_time::Entry{T} = Entry(T, "year", "Duration of construction"; default=8.0)
    fixed_charge_rate::Entry{T} = Entry(T, "-", "Constant dollar fixed charge rate"; default=0.078)
    initial_cost_blanket::Entry{T} = Entry(T, "\$M", "Cost of initial blanket"; default=200.0)
    initial_cost_divertor::Entry{T} = Entry(T, "\$M", "Cost of initial divertor"; default=8.0)
    divertor_fluence_lifetime::Entry{T} = Entry(T, "MW*yr/m^2", "Divertor fluence over its lifetime"; default=10.0)
    blanket_fluence_lifetime::Entry{T} = Entry(T, "MW*yr/m^2", "Blanket fluence over its lifetime"; default=15.0)
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
    year_assessed::Union{Missing,Int}
    year::Union{Missing,Int}
end

function DollarAdjust(par::FUSEparameters__ActorCosting)
    return DollarAdjust(par.future_inflation_rate, par.construction_start_year, missing, missing)
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

    # from old dollars to current dollars
    if CPI.Year[1] <= da.year_assessed <= CPI.Year[end-1]
        CPI_past_year = CPI[findfirst(x -> x==da.year_assessed, CPI.Year), "Year Avg"]
        val_today = CPI[end, "Year Avg"] ./ CPI_past_year .* dollars
    elseif da.year_assessed == CPI.Year[end]
        val_today = dollars
    else
        @warn "Inflation data not available for $(da.year_assessed)"
        val_today = dollars
    end

    # inflate to the year of start of construction
    if da.construction_start_year < CPI.Year[1]
        error("Cannot translate cost earlier than $(CPI.Year[1])")
    elseif da.construction_start_year <= CPI.Year[end]
        CPI_const_year = CPI[findfirst(x -> x==da.construction_start_year, CPI.Year), "Year Avg"]
        value = CPI_const_year ./ CPI[end, "Year Avg"] .* val_today
    else
        n_years = da.construction_start_year - CPI.Year[end]
        value = val_today * ((1.0 + da.future_inflation_rate)^n_years)
    end

    # wipe out year_assessed each time to force developer to enter of `da.year_assessed` and avoid using the wrong year 
    da.year_assessed = missing

    return value
end

"""
    ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the cost of building, operating, and recommission the fusion power plant.

!!! note 
    Stores data in `dd.costing`
"""
function ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCosting(dd, act.ActorCosting; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorCosting)
    par = actor.par
    dd = actor.dd
    cst = dd.costing

    empty!(cst)

    cst.construction_start_year = par.construction_start_year
    cst.plant_lifetime = par.plant_lifetime

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
