#= ================= =#
#  ActorCostingARIES  #
#= ================= =#

Base.@kwdef mutable struct FUSEparameters__ActorCostingARIES{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    land_space::Entry{T} = Entry{T}("acres", "Plant site space required"; default=1000.0)
    building_volume::Entry{T} = Entry{T}("m^3", "Volume of the tokmak building"; default=140.0e3)
    interest_rate::Entry{T} = Entry{T}("-", "Annual interest rate fraction of direct capital cost"; default=0.05)
    indirect_cost_rate::Entry{T} =
        Entry{T}("-", "Indirect cost associated with construction, equipment, services, engineering construction management and owners cost"; default=0.4)
    escalation_fraction::Entry{T} = Entry{T}("-", "Yearly escalation fraction based on risk assessment"; default=0.05)
    blanket_lifetime::Entry{T} = Entry{T}("year", "Lifetime of the blanket"; default=6.8)
end

mutable struct ActorCostingARIES{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCostingARIES{P}
    function ActorCostingARIES(dd::IMAS.dd{D}, par::FUSEparameters__ActorCostingARIES{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorCostingARIES)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorCostingARIES(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates costing based on ARIES cost account documentation https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf

!!! note

    Stores data in `dd.costing`
"""
function ActorCostingARIES(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCostingARIES(kw...)
    actor = ActorCostingARIES(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorCostingARIES)
    dd = actor.dd
    par = actor.par
    cst = dd.costing

    cost_direct = cst.cost_direct_capital
    cost_ops = cst.cost_operations
    cost_decom = cst.cost_decommissioning

    da = DollarAdjust(dd)

    ###### Direct Capital ######

    ### Tokamak

    # build layers
    tokamak = resize!(cost_direct.system, "name" => "tokamak")
    for layer in dd.build.layer
        if layer.side == Int(_lfs_)
            continue # avoid double counting of hfs and lfs layers
        elseif layer.type == Int(_oh_)
            continue # avoid double counting of oh
        end
        c = cost_direct_capital_ARIES(layer, cst, da)
        if c > 0.0
            sub = resize!(tokamak.subsystem, "name" => replace(layer.name, r"^hfs " => ""))
            sub.cost = c
        end
    end

    # PF coils
    for (name, c) in cost_direct_capital_ARIES(dd.pf_active, da, cst)
        sub = resize!(tokamak.subsystem, "name" => name)
        sub.cost = c
    end

    # Heating and current drive
    for hcd in vcat(dd.ec_launchers.beam, dd.ic_antennas.antenna, dd.lh_antennas.antenna, dd.nbi.unit)
        c = cost_direct_capital_ARIES(hcd, da)
        if c > 0.0
            sub = resize!(tokamak.subsystem, "name" => uppercase(hcd.name))
            sub.cost = c
        end
    end

    ### Facility

    sys = resize!(cost_direct.system, "name" => "facility")

    power_thermal, power_electric_generated, power_electric_net = bop_powers(dd.balance_of_plant)

    for item in vcat(:land, :buildings, :hot_cell, :heat_transfer_loop_materials, :balance_of_plant_equipment, :fuel_cycle_rad_handling)
        sub = resize!(sys.subsystem, "name" => replace(string(item), "_" => " "))
        if item == :land
            sub.cost = cost_direct_capital_ARIES(item, par.land_space, power_electric_generated, da)
        elseif item == :buildings
            sub.cost = cost_direct_capital_ARIES(item, par.building_volume,
                par.land_space, power_electric_generated,
                power_thermal, power_electric_net, da)
        elseif item == :hot_cell
            sub.cost = cost_direct_capital_ARIES(item, par.building_volume, da)
        elseif item == :heat_transfer_loop_materials
            sub.cost = cost_direct_capital_ARIES(item, power_thermal, da)
        elseif item == :balance_of_plant_equipment
            sub.cost = cost_direct_capital_ARIES(item, power_thermal, power_electric_generated, da, dd)
        elseif item == :fuel_cycle_rad_handling
            sub.cost = cost_direct_capital_ARIES(item, power_thermal, power_electric_net, da)
        else
            sub.cost = cost_direct_capital_ARIES(item, da)
        end
    end

    ###### Operations ######
    sys = resize!(cost_ops.system, "name" => "tritium handling")
    sys.yearly_cost = cost_operations_ARIES(:tritium_handling, da)

    sys = resize!(cost_ops.system, "name" => "maintenance and operators")
    sys.yearly_cost = cost_operations_ARIES(:operation_maintenance, power_electric_generated, da)

    sys = resize!(cost_ops.system, "name" => "replacements")
    for item in (:blanket_replacement,)
        sub = resize!(sys.subsystem, "name" => string(item))
        if item == :blanket_replacement
            blanket_cost = findfirst(item -> item.name == "blanket", tokamak.subsystem)
            sub.yearly_cost = cost_operations_ARIES(:blanket_replacement, blanket_cost, par.blanket_lifetime, da)
        else
            sub.yearly_cost = cost_operations_ARIES(item, da)
        end
    end

    ###### Decomissioning ######
    sys = resize!(cost_decom.system, "name" => "decommissioning")
    sys.cost = cost_decomissioning_ARIES(:decom_wild_guess, cst.plant_lifetime, da)

    ###### Levelized Cost Of Electricity 
    capital_cost_rate = par.interest_rate / (1.0 - (1.0 + par.interest_rate)^(-1.0 * cst.plant_lifetime))
    cst.cost_lifetime = 0.0
    for year in 1:cst.plant_lifetime
        yearly_cost = (capital_cost_rate * cost_direct.cost + cost_ops.yearly_cost + cost_decom.cost / cst.plant_lifetime)
        cst.cost_lifetime += (1.0 + par.escalation_fraction) * (1.0 + par.indirect_cost_rate) * yearly_cost
    end
    cst.levelized_CoE = (cst.cost_lifetime * 1E6) / (cst.plant_lifetime * 24 * 365 * power_electric_net / 1e3 * cst.availability)

    return actor
end

function _finalize(actor::ActorCostingARIES)
    # sort system/subsystem by their costs
    sort!(actor.dd.costing.cost_direct_capital.system; by=x -> x.cost, rev=true)
    for sys in actor.dd.costing.cost_direct_capital.system
        sort!(sys.subsystem; by=x -> x.cost, rev=true)
    end

    return actor
end

#= =================== =#
#  direct capital cost  #
#= =================== =#
"""
    cost_direct_capital_ARIES(layer::IMAS.build__layer, cst::IMAS.costing, da::DollarAdjust)

Capital cost for each layer in the build
"""
function cost_direct_capital_ARIES(layer::IMAS.build__layer, cst::IMAS.costing, da::DollarAdjust)
    da.year_assessed = 2016
    if layer.type == Int(_oh_)
        return 0.0
    elseif layer.type == Int(_tf_)
        build = IMAS.parent(IMAS.parent(layer))
        cost = layer.volume * unit_cost(build.tf.technology, cst)
        return future_dollars(cost, da)
    elseif layer.type == Int(_shield_)
        cost = layer.volume * 0.29  # $M/m^3
        return future_dollars(cost, da)
    elseif layer.type == Int(_blanket_)
        cost = layer.volume * 0.75  # $M/m^3
        return future_dollars(cost, da)
    elseif layer.type ∈ (Int(_wall_), Int(_vessel_), Int(_cryostat_))
        cost = layer.volume * 0.36  # $M/m^3
        return future_dollars(cost, da)
    else
        cost = layer.volume * unit_cost(Material(layer.material), cst)
        return future_dollars(cost, da)
    end
end

"""
    cost_direct_capital_ARIES(ecl::IMAS.ec_launchers__beam, da::DollarAdjust)

Capital cost for each EC launcher is proportional to its power

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(ecl::IMAS.ec_launchers__beam, da::DollarAdjust)
    da.year_assessed = 2009
    unit_cost = 3.0 # $/W
    cost = ecl.available_launch_power / 1E6 * unit_cost
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(ica::IMAS.ic_antennas__antenna, da::DollarAdjust)

Capital cost for each IC antenna is proportional to its power

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(ica::IMAS.ic_antennas__antenna, da::DollarAdjust)
    da.year_assessed = 2009   #Table 8 in ARIES
    unit_cost = 1.64 # $/W
    cost = ica.available_launch_power / 1E6 * unit_cost
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(lha::IMAS.lh_antennas__antenna, da::DollarAdjust)

Capital cost for each LH antenna is proportional to its power

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(lha::IMAS.lh_antennas__antenna, da::DollarAdjust)
    da.year_assessed = 2009
    unit_cost = 2.13 # $/W
    cost = lha.available_launch_power / 1E6 * unit_cost
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(nbu::IMAS.nbi__unit, da::DollarAdjust)

Capital cost for each NBI unit is proportional to its power

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(nbu::IMAS.nbi__unit, da::DollarAdjust)
    da.year_assessed = 2009
    unit_cost = 4.93 # $/W
    cost = nbu.available_launch_power / 1E6 * unit_cost
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(pf_active::IMAS.pf_active, da::DollarAdjust, cst::IMAS.costing)
"""
function cost_direct_capital_ARIES(pf_active::IMAS.pf_active, da::DollarAdjust, cst::IMAS.costing)
    dd = IMAS.top_dd(pf_active)
    c = Dict("OH" => 0.0, "PF" => 0.0)
    for coil in pf_active.coil
        if IMAS.is_ohmic_coil(coil)
            c["OH"] += cost_direct_capital_ARIES(coil, dd.build.oh.technology, da, cst)
        else
            c["PF"] += cost_direct_capital_ARIES(coil, dd.build.pf_active.technology, da, cst)
        end
    end
    return c
end

"""
    cost_direct_capital_ARIES(coil::IMAS.pf_active__coil, technology::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology}, da::DollarAdjust, cst::IMAS.costing)
"""
function cost_direct_capital_ARIES(
    coil::IMAS.pf_active__coil,
    technology::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology},
    da::DollarAdjust,
    cst::IMAS.costing
)
    da.year_assessed = 2016
    cost = IMAS.volume(coil) * unit_cost(technology, cst)
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(::Type{Val{:land}}, land::Real, power_electric_net::Real, da::DollarAdjust)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(::Type{Val{:land}}, land::Real, power_electric_net::Real, da::DollarAdjust)
    da.year_assessed = 2009 #pg. 8, ARIES report
    power_electric_net = power_electric_net / 1E6
    if power_electric_net > 0.0
        cost = 1.2 * land * 20.0e-3 * (power_electric_net / 1000.0)^0.3
    else
        cost = 1.2 * land * 20.0e-3
    end
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(::Type{Val{:buildings}}, land::Real, building_volume::Real, power_electric_generated::Real, power_thermal::Real, power_electric_net::Real, da::DollarAdjust)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(
    ::Type{Val{:buildings}},
    land::Real,
    building_volume::Real,
    power_electric_generated::Real,
    power_thermal::Real,
    power_electric_net::Real,
    da::DollarAdjust
)
    da.year_assessed = 2009
    power_electric_generated = power_electric_generated / 1E6
    power_thermal = power_thermal / 1E6
    power_electric_net = power_electric_net / 1E6
    cost = 27.0 * (land / 1000.0)^0.2 # site
    cost += 111.661 * (building_volume / 80.0e3)^0.62 # tokamak building
    cost += 78.9 * (power_electric_generated / 1246)^0.5 # turbine/generator
    cost += 16.804 * ((power_thermal - power_electric_net) / 1860.0)^0.5 # heat rejection
    cost += 22.878 * (power_electric_net / 1000.0)^0.3 # electrical equipment
    cost += 21.96 * (power_electric_generated / 1246)^0.3  # Plant auxiliary systems
    cost += 4.309 * (power_electric_net / 1000.0)^0.3 # power core service building
    cost += 1.513 * (power_electric_net / 1000.0)^0.3  # service water
    cost += 25.0 * (power_thermal / 1759.0)^0.3  # tritium plant and systems
    cost += 7.11 # control room
    cost += 4.7 * (power_electric_net / 1000.0)^0.3  # On-site AC power supply building
    cost += 2.0 # administrative
    cost += 2.0 # site service
    cost += 2.09 # cyrogenic and inert gas storage
    cost += 0.71 # security
    cost += 4.15  # ventilation stack
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(::Type{Val{:hot_cell}}, building_volume, da::DollarAdjust)

NOTE: https://www.iter.org/mach/HotCell
"""
function cost_direct_capital_ARIES(::Type{Val{:hot_cell}}, building_volume::Real, da::DollarAdjust)
    da.year_assessed = 2009
    cost = 0.34 * 111.661 * (building_volume / 80.0e3)^0.62
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(::Type{Val{:heat_transfer_loop_materials}}, power_thermal::Real, da::DollarAdjust)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf (warning uses LiPb for blanket)
"""
function cost_direct_capital_ARIES(::Type{Val{:heat_transfer_loop_materials}}, power_thermal::Real, da::DollarAdjust)
    da.year_assessed = 2009  # Table 15 ARIES report
    power_thermal = power_thermal / 1E6
    cost = 50.0 * (power_thermal / 2000.0)^0.55 # water
    cost += 125.0 * (power_thermal / 2000.0)^0.55 # LiPb
    cost += 110.0 * (power_thermal / 2000.0)^0.55 # He
    cost += 0.01 * power_thermal # NbIHX
    cost += 50.0 * (power_thermal / 2000.0)^0.55 # Na
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(::Type{Val{:balance_of_plant_equipment}}, power_thermal::Real, power_electric_generated::Real, da::DollarAdjust, dd::IMAS.dd)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(::Type{Val{:balance_of_plant_equipment}}, power_thermal::Real, power_electric_generated::Real, da::DollarAdjust, dd::IMAS.dd)
    da.year_assessed = 2009
    power_thermal = power_thermal / 1E6
    power_electric_generated = power_electric_generated / 1E6
    bop = dd.balance_of_plant

    if contains(lowercase(bop.power_plant.power_cycle_type), "rankine")
        cost = 350.0 * (power_thermal / 2620.0)^0.7 # Turbine equipment
    elseif contains(lowercase(bop.power_plant.power_cycle_type), "brayton")
        cost = 360.0 * (power_thermal / 2000.0)^0.8
        if !isnan(bop.thermal_efficiency_cycle[1])
            cost *= bop.thermal_efficiency_cycle[1] / 0.6
        end
    end
    cost += 182.98 * (power_electric_generated / 1200.0)^0.5 # Electrical plant equipment
    cost += 87.52 * ((power_thermal - power_electric_generated) / 2300.0) # Heat rejection equipment
    cost += 88.89 * (power_electric_generated / 1200.0)^0.6 # Miscellaneous equipment
    return future_dollars(cost, da)
end

"""
    cost_direct_capital_ARIES(::Type{Val{:fuel_cycle_rad_handling}}, power_thermal::Real, power_electric_net::Real, da::DollarAdjust)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf (warning uses LiPb for blanket)
"""
function cost_direct_capital_ARIES(::Type{Val{:fuel_cycle_rad_handling}}, power_thermal::Real, power_electric_net::Real, da::DollarAdjust)
    da.year_assessed = 2009
    power_thermal = power_thermal / 1E6
    power_electric_net = power_electric_net / 1E6
    cost = 15.0 * (power_thermal / 1758.0)^0.85 # radioactive material treatment and management
    cost += 70.0 * (power_thermal / 1758.0)^0.8 # Fuel handling and storage
    cost += 100.0 * (power_electric_net / 2000.0)^0.55 # Hot cell maintenance
    cost += 60.0 # Instrumentation and Control
    cost += 8.0 * (power_thermal / 1000.0)^0.8 # Misc power core equipment
    return future_dollars(cost, da)
end

#= ====================== =#
#  yearly operations cost  #
#= ====================== =#
"""
    cost_operations_ARIES(::Type{Val{:operation_maintenance}}, power_electric_generated::Real, da::DollarAdjust)

Yearly cost for maintenance [\$M/year]

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_operations_ARIES(::Type{Val{:operation_maintenance}}, power_electric_generated::Real, da::DollarAdjust)
    da.year_assessed = 2009
    power_electric_generated = power_electric_generated / 1E6

    if power_electric_generated == 0.0
        cost = 80.0 # estimate based on Table 36 in ARIES
    else
        cost = 80.0 * (power_electric_generated / 1200.0)^0.5
    end
    return future_dollars(cost, da)
end

"""
    cost_operations_ARIES(::Type{Val{:tritium_handling}}, da::DollarAdjust)

Yearly cost for tritium_handling [\$M/year]

!!!!WRONG!!!! Needs estiamte
"""
function cost_operations_ARIES(::Type{Val{:tritium_handling}}, da::DollarAdjust)
    da.year_assessed = 2013
    cost = 1.0
    return future_dollars(cost, da)
end

"""
    cost_operations_ARIES(::Type{Val{:blanket_replacement}}, cost_blanket::Real, blanket_lifetime::Real, da::DollarAdjust)

Yearly cost for blanket replacement [\$M/year]
"""
function cost_operations_ARIES(::Type{Val{:blanket_replacement}}, cost_blanket::Union{Real,Nothing}, blanket_lifetime::Real, da::DollarAdjust)
    if cost_blanket === nothing
        cost = 0.0
    else
        da.year_assessed = Dates.year(Dates.now()) # Assume that the user will give you the cost of the blanket in the dollars of their current year 
        cost = cost_blanket / blanket_lifetime
    end
    return future_dollars(cost, da)
end

#= =================== =#
#  Decomissioning cost  #
#= =================== =#
"""
    cost_decomissioning_ARIES(::Type{Val{:decom_wild_guess}}, plant_lifetime::Real, da::DollarAdjust)

Cost to decommission the plant [\$M]

LIKELY NEEDS FIXING
"""
function cost_decomissioning_ARIES(::Type{Val{:decom_wild_guess}}, plant_lifetime::Real, da::DollarAdjust)
    da.year_assessed = 2009  # pg. 94 ARIES
    unit_cost = 2.76 # [$M/year] from GASC
    cost = unit_cost * plant_lifetime
    return future_dollars(cost, da)
end
