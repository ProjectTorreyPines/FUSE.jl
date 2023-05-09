using CSV
using DataFrames
import Memoize
import Dates

#= ================================= =#
#  Learning rate for HTS - from GASC  #
#= ================================= =#
#assumes that production volume of HTS has increased by a factor of 10 and that there is a 15% decrease in cost for every doubling of production volume
function cost_multiplier(production_increase::Real, learning_rate::Real)
    return production_increase^(log(learning_rate)/log(2))
end

#= ============== =#
#  materials cost  #
#= ============== =#
#NOTE: material should be priced by Kg
#NOTE: if something is priced by m^3 then it is for a specific part already
function unit_cost(material::AbstractString, production_increase::Real, learning_rate::Real)
    if material == "Vacuum"
        return 0.0 # $M/m^3
    elseif material == "ReBCO"
        return (87.5 / 2) * cost_multiplier(production_increase, learning_rate) # $M/m^3
    elseif material == "Nb3Sn"
        return 1.66 # $M/m^3
    elseif contains(lowercase(material), "steel")
        return 0.36 # $M/m^3
    elseif material == "Tungsten"
        return 0.36 # $M/m^3
    elseif material == "Copper"
        return 0.5 # $M/m^3
    elseif material == "Water, Liquid"
        return 0.0 # $M/m^3
    elseif material == "lithium-lead"
        return 0.75 # $M/m^3
    elseif material == "FLiBe"
        return 0.75 * 3 # $M/m^3
    elseif contains(lowercase(material), "plasma")
        return 0.0 # $M/m^3
    else
        error("Material `$material` has no price \$M/m³")
    end
end

#= ====================== =#
#  materials cost - coils  #
#= ====================== =#
function unit_cost(coil_tech::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology}, production_increase::Real, learning_rate::Real)
    if coil_tech.material == "Copper"
        return unit_cost("Copper", production_increase, learning_rate)
    else
        fraction_cable = 1 - coil_tech.fraction_stainless - coil_tech.fraction_void
        fraction_SC = fraction_cable * coil_tech.ratio_SC_to_copper / (1 + coil_tech.ratio_SC_to_copper)
        fraction_copper = fraction_cable - fraction_SC
        return (coil_tech.fraction_stainless * unit_cost("Steel, Stainless 316", production_increase, learning_rate) + fraction_copper * unit_cost("Copper", production_increase, learning_rate) + fraction_SC * unit_cost(coil_tech.material, production_increase, learning_rate))
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

function DollarAdjust(dd::IMAS.dd)
    return DollarAdjust(dd.costing.future_inflation_rate, dd.costing.construction_start_year, missing, missing)
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

#= ================== =#
#  Dispatch on symbol  #
#= ================== =#

#Sheffield

function cost_direct_capital_Sheffield(item::Symbol, args...; kw...)
    return cost_direct_capital_Sheffield(Val{item}, args...; kw...)
end

function cost_ops_Sheffield(item::Symbol, args...; kw...)
    return cost_ops_Sheffield(Val{item}, args...; kw...)
end

function cost_fuel_Sheffield(item::Symbol, args...; kw...)
    return cost_fuel_Sheffield(Val{item}, args...; kw...)
end

function cost_operations_Sheffield(item::Symbol, args...; kw...)
    return cost_operations_Sheffield(Val{item}, args...; kw...)
end

function cost_decomissioning_Sheffield(item::Symbol, args...; kw...)
    return cost_decomissioning_Sheffield(Val{item}, args...; kw...)
end

#ARIES

function cost_direct_capital_ARIES(item::Symbol, args...; kw...)
    return cost_direct_capital_ARIES(Val{item}, args...; kw...)
end

function cost_operations_ARIES(item::Symbol, args...; kw...)
    return cost_operations_ARIES(Val{item}, args...; kw...)
end

function cost_decomissioning_ARIES(item::Symbol, args...; kw...)
    return cost_decomissioning_ARIES(Val{item}, args...; kw...)
end