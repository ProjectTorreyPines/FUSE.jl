import CSV
import DataFrames
import Memoize
import Dates
import Measurements: value, uncertainty
import FusionMaterials
import FusionMaterials: Material

#= ================================= =#
#  Learning rate for HTS - from GASC  #
#= ================================= =#
function cost_multiplier(production_increase::Real, learning_rate::Real)
    return production_increase^(log(learning_rate) / log(2))
end

#= ============== =#
#  materials cost  #
#= ============== =#
function unit_cost(material::Material, cst::IMAS.costing)
    cost_per_unit_volume = material.cost_m3 * material.manufacturing_multiplier / 1e6 # costs in $M/m^3

    if material.name == "rebco"
        production_increase = cst.future.learning.hts.production_increase
        learning_rate = cst.future.learning.hts.learning_rate
        cost_per_unit_volume = cost_per_unit_volume * cost_multiplier(production_increase, learning_rate)
    end

    return cost_per_unit_volume
end

#= ====================== =#
#  materials cost - coils  #
#= ====================== =#
function unit_cost(coil_tech::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology}, cst::IMAS.costing)
    if coil_tech.material == "copper"
        return unit_cost(Material(:copper), cst)
    else
        fraction_cable = 1.0 - coil_tech.fraction_steel - coil_tech.fraction_void
        fraction_SC = fraction_cable * coil_tech.ratio_SC_to_copper / (1 + coil_tech.ratio_SC_to_copper)
        fraction_copper = fraction_cable - fraction_SC
        cost =
            coil_tech.fraction_steel * unit_cost(Material(:steel), cst) +
            fraction_copper * unit_cost(Material(:copper), cst) +
            fraction_SC * unit_cost(Material(coil_tech.material), cst)
        return cost
    end
end

#= ==================== =#
#  Inflation Adjustment  #
#= ==================== =#
mutable struct DollarAdjust{T<:Real}
    future_inflation_rate::T
    construction_start_year::T
    year_assessed::Union{Missing,Int}
    year::Union{Missing,Int}
end

function DollarAdjust(dd::IMAS.dd{T}) where {T<:Real}
    return DollarAdjust{T}(dd.costing.future.inflation_rate, dd.costing.construction_start_year, missing, missing)
end
Base.eltype(::Type{DollarAdjust{T}}) where {T} = T

"""
    load_inflation_rate()

Inflation of past costs according to U.S. Bureau of Labor Statistics CPI data

Source: https://data.bls.gov/cgi-bin/surveymost?cu (Dataset CUUR0000SA0)

NOTE: To update simply copy-paste the table at the website in the "CPI.csv" file
"""
Memoize.@memoize function load_inflation_rate()
    months = [:Jan, :Feb, :Mar, :Apr, :May, :Jun, :Jul, :Aug, :Sep, :Oct, :Nov, :Dec]
    types = Dict(
        :Year => Int64,
        :HALF1 => Union{Float64,Missing},
        :HALF2 => Union{Float64,Missing}
    )
    for month in months
        types[month] = Union{Float64,Missing}
    end

    csv_loc = abspath(joinpath(@__DIR__, "CPI.csv"))

    CPI = DataFrames.DataFrame(CSV.File(csv_loc; types, missingstring=" ", silencewarnings=true))
    DataFrames.select!(CPI, Not(:HALF1))
    DataFrames.select!(CPI, Not(:HALF2))

    CPI[!, "Year Avg"] = [sum([x for x in row[months] if x !== missing]) / length([x for x in row[months] if x !== missing]) for row in eachrow(CPI)]

    return CPI
end

"""
    future_dollars(dollars::Real, da::DollarAdjust)

Adjusts costs assessed in a previous year to current or future dollar amount
"""
function future_dollars(dollars::Real, da::DollarAdjust{T})::T where {T<:Real} 

    if dollars == 0.0
        return zero(T)
    end

    CPI = load_inflation_rate()

    # from old dollars to current dollars
    if CPI.Year[1] <= da.year_assessed <= CPI.Year[end]
        CPI_past_year = CPI[findfirst(x -> x == da.year_assessed, CPI.Year), "Year Avg"]
        val_today = CPI[end, "Year Avg"] ./ CPI_past_year .* dollars
    elseif da.year_assessed == CPI.Year[end] + 1
        # allow one year of slack before raising warning that inflation data is not available
        val_today = dollars
    else
        @warn "Inflation data not available for $(da.year_assessed)"
        val_today = dollars
    end

    # inflate to the year of start of construction
    cpi_year_min = CPI.Year[1]
    cpi_year_max = CPI.Year[end]
    construction_start_year = da.construction_start_year
    construction_year_min, construction_year_max = year_bounds(construction_start_year)

    if construction_year_min < cpi_year_min
        error("Cannot translate construction_start_year=$(construction_start_year) earlier than $(cpi_year_min)")
    elseif construction_year_min <= cpi_year_max && construction_year_max > cpi_year_max
        error("construction_start_year=$(construction_start_year) overlaps CPI range end $(cpi_year_max); use a fully historical or fully future year range")
    elseif construction_year_max <= cpi_year_max
        cpi_construction_year = cpi_year_average(CPI, construction_start_year)
        value = cpi_construction_year ./ CPI[end, "Year Avg"] .* val_today
    else
        n_years = construction_start_year - cpi_year_max
        value = val_today * ((1.0 + da.future_inflation_rate)^n_years)
    end

    # wipe out year_assessed each time to force developer to enter of `da.year_assessed` and avoid using the wrong year 
    da.year_assessed = missing

    return T(value)
end

function year_bounds(year::Real)
    year_real = float(year)
    return year_real, year_real
end

function year_bounds(year::Measurement)
    y0 = value(year)
    σ = uncertainty(year)
    return y0 - σ, y0 + σ
end

function cpi_year_average(CPI::DataFrames.DataFrame, year::Real)
    cpi_year_min = CPI.Year[1]
    cpi_year_max = CPI.Year[end]
    if !(cpi_year_min <= year <= cpi_year_max)
        error("Requested year=$(year) is outside CPI range [$cpi_year_min, $cpi_year_max]")
    end

    year_lo = floor(Int, year)
    year_hi = ceil(Int, year)
    idx_lo = findfirst(==(year_lo), CPI.Year)
    idx_hi = findfirst(==(year_hi), CPI.Year)
    @assert idx_lo !== nothing && idx_hi !== nothing "CPI table is missing expected years for interpolation"

    cpi_lo = CPI[idx_lo, "Year Avg"]
    cpi_hi = CPI[idx_hi, "Year Avg"]

    if year_lo == year_hi
        return cpi_lo
    end

    weight = year - year_lo
    return (1.0 - weight) * cpi_lo + weight * cpi_hi
end

function cpi_year_average(CPI::DataFrames.DataFrame, year::Measurement)
    cpi_year_min = CPI.Year[1]
    cpi_year_max = CPI.Year[end]
    year_min, year_max = year_bounds(year)
    if !(cpi_year_min <= year_min <= year_max <= cpi_year_max)
        error("Requested year=$(year) is outside CPI range [$cpi_year_min, $cpi_year_max]")
    end

    year_nominal = value(year)
    year_uncertainty = uncertainty(year)
    cpi_nominal = cpi_year_average(CPI, year_nominal)
    cpi_gradient = cpi_year_average_gradient(CPI, year_nominal)
    cpi_uncertainty = abs(cpi_gradient) * year_uncertainty

    return cpi_nominal ± cpi_uncertainty
end

function cpi_year_average_gradient(CPI::DataFrames.DataFrame, year::Real)
    cpi_year_min = CPI.Year[1]
    cpi_year_max = CPI.Year[end]
    if !(cpi_year_min <= year <= cpi_year_max)
        error("Requested year=$(year) is outside CPI range [$cpi_year_min, $cpi_year_max]")
    end

    # Piecewise-linear CPI interpolation slope d(CPI)/d(year).
    # At exact integer years (kinks), use the right derivative except at the last year.
    year_lo = floor(Int, year)
    year_hi = ceil(Int, year)
    if year_lo == year_hi
        if year_lo == cpi_year_max
            year_lo = cpi_year_max - 1
            year_hi = cpi_year_max
        else
            year_hi = year_lo + 1
        end
    end

    idx_lo = findfirst(==(year_lo), CPI.Year)
    idx_hi = findfirst(==(year_hi), CPI.Year)
    @assert idx_lo !== nothing && idx_hi !== nothing "CPI table is missing expected years for gradient evaluation"

    cpi_lo = CPI[idx_lo, "Year Avg"]
    cpi_hi = CPI[idx_hi, "Year Avg"]

    return (cpi_hi - cpi_lo) / (year_hi - year_lo)
end

#= ================== =#
#  Dispatch on symbol  #
#= ================== =#

#Sheffield (and GASC)

function cost_direct_capital_Sheffield(item::Symbol, args...; kw...)
    return cost_direct_capital_Sheffield(Val(item), args...; kw...)
end

function cost_ops_Sheffield(item::Symbol, args...; kw...)
    return cost_ops_Sheffield(Val(item), args...; kw...)
end

function cost_fuel_Sheffield(item::Symbol, args...; kw...)
    return cost_fuel_Sheffield(Val(item), args...; kw...)
end

function cost_operations_Sheffield(item::Symbol, args...; kw...)
    return cost_operations_Sheffield(Val(item), args...; kw...)
end

function cost_decomissioning_Sheffield(item::Symbol, args...; kw...)
    return cost_decomissioning_Sheffield(Val(item), args...; kw...)
end

#ARIES

function cost_direct_capital_ARIES(item::Symbol, args...; kw...)
    return cost_direct_capital_ARIES(Val(item), args...; kw...)
end

function cost_operations_ARIES(item::Symbol, args...; kw...)
    return cost_operations_ARIES(Val(item), args...; kw...)
end

function cost_decomissioning_ARIES(item::Symbol, args...; kw...)
    return cost_decomissioning_ARIES(Val(item), args...; kw...)
end

#= ========== =#
#  BOP powers  #
#= ========== =#
"""
    bop_powers(bop::IMAS.balance_of_plant)

Returns maximum total_useful_heat_power, power_electric_generated, power_electric_net over time, setting them to 0.0 if negative, NaN or missing
"""
function bop_powers(bop::IMAS.balance_of_plant)
    if ismissing(bop.power_plant, :total_heat_supplied)
        total_useful_heat_power = 0.0
    else
        total_useful_heat_power = max(0.0, maximum(x -> isnan(x) ? -Inf : x, bop.power_plant.total_heat_supplied))
    end

    if ismissing(bop.power_plant, :power_electric_generated)
        power_electric_generated = 0.0
    else
        power_electric_generated = max(0.0, maximum(x -> isnan(x) ? -Inf : x, bop.power_plant.power_electric_generated))
    end

    if ismissing(bop, :power_electric_net)
        power_electric_net = 0.0
    else
        power_electric_net = max(0.0, maximum(x -> isnan(x) ? -Inf : x, bop.power_electric_net))
    end

    return (total_useful_heat_power=total_useful_heat_power, power_electric_generated=power_electric_generated, power_electric_net=power_electric_net)
end
