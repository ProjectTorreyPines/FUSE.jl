#= ================== =#
#  Dispatch on symbol  #
#= ================== =#
function cost_direct_capital_ARIES(item::Symbol, args...; kw...)
    return cost_direct_capital_ARIES(Val{item}, args...; kw...)
end

function cost_operations_ARIES(item::Symbol, args...; kw...)
    return cost_operations_ARIES(Val{item}, args...; kw...)
end

function cost_decomissioning_ARIES(item::Symbol, args...; kw...)
    return cost_decomissioning_ARIES(Val{item}, args...; kw...)
end

#= ============== =#
#  materials cost  #
#= ============== =#

function unit_cost(coil_tech::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology})
    if coil_tech.material == "Copper"
        return unit_cost("Copper")
    else
        fraction_cable = 1 - coil_tech.fraction_stainless - coil_tech.fraction_void
        fraction_SC = fraction_cable * coil_tech.ratio_SC_to_copper
        fraction_copper = fraction_cable - fraction_SC
        return (coil_tech.fraction_stainless * unit_cost("Steel, Stainless 316") + fraction_copper * unit_cost("Copper") + fraction_SC * unit_cost(coil_tech.material))
    end
end

#= =================== =#
#  direct capital cost  #
#= =================== =#
"""
    cost_direct_capital_ARIES(layer::IMAS.build__layer)

Capital cost for each layer in the build
"""
function cost_direct_capital_ARIES(layer::IMAS.build__layer)
    if layer.type == Int(_oh_)
        build = IMAS.parent(IMAS.parent(layer))
        return layer.volume * unit_cost(build.oh.technology)
    elseif layer.type == Int(_tf_)
        build = IMAS.parent(IMAS.parent(layer))
        return layer.volume * unit_cost(build.tf.technology)
    elseif layer.type == Int(_shield_)
        return layer.volume * 0.29  # $M/m^3
    elseif layer.type == Int(_blanket_)
        return layer.volume * 0.75  # $M/m^3
    elseif layer.type ∈ [Int(_wall_), Int(_vessel_), Int(_cryostat_)]
        return layer.volume * 0.36  # $M/m^3
    else
        return layer.volume * unit_cost(layer.material)
    end
end

"""
    cost_direct_capital_ARIES(ecl::IMAS.ec_launchers__beam)

Capital cost for each EC launcher is proportional to its power

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(ecl::IMAS.ec_launchers__beam)
    unit_cost = 3.0 # $/W
    return ecl.available_launch_power / 1E6 * unit_cost
end

"""
    cost_direct_capital_ARIES(ica::IMAS.ic_antennas__antenna)

Capital cost for each IC antenna is proportional to its power

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(ica::IMAS.ic_antennas__antenna)
    unit_cost = 1.64 # $/W
    return ica.available_launch_power / 1E6 * unit_cost
end

"""
    cost_direct_capital_ARIES(lha::IMAS.lh_antennas__antenna)

Capital cost for each LH antenna is proportional to its power

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(lha::IMAS.lh_antennas__antenna)
    unit_cost = 2.13 # $/W
    return lha.available_launch_power / 1E6 * unit_cost
end

"""
    cost_direct_capital_ARIES(nbu::IMAS.nbi__unit)

Capital cost for each NBI unit is proportional to its power

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(nbu::IMAS.nbi__unit)
    unit_cost = 4.93 # $/W
    return nbu.available_launch_power / 1E6 * unit_cost
end

"""
    cost_direct_capital_ARIES(pf_active::IMAS.pf_active)
"""
function cost_direct_capital_ARIES(pf_active::IMAS.pf_active)
    dd = IMAS.top_dd(pf_active)
    c = Dict("OH" => 0.0, "PF" => 0.0)
    for coil in pf_active.coil
        if coil.name == "OH"
            c["OH"] += cost_direct_capital_ARIES(coil, dd.build.oh.technology)
        else
            c["PF"] += cost_direct_capital_ARIES(coil, dd.build.pf_active.technology)
        end
    end
    return c
end

"""
    cost_direct_capital_ARIES(coil::IMAS.pf_active__coil, technology::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology})

"""
function cost_direct_capital_ARIES(coil::IMAS.pf_active__coil, technology::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology})
    return IMAS.volume(coil) * unit_cost(technology)
end

"""
    cost_direct_capital_ARIES(::Type{Val{:land}}, land::Real, power_electric_net::Real)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(::Type{Val{:land}}, land::Real, power_electric_net::Real)
    power_electric_net = power_electric_net / 1E6
    return 1.2 * land * 20.0e-3 * (power_electric_net / 1000.0)^0.3
end

"""
    cost_direct_capital_ARIES(::Type{Val{:buildings}}, land::Real, building_volume::Real, power_electric_generated::Real, power_thermal::Real, power_electric_net::Real)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(::Type{Val{:buildings}}, land::Real, building_volume::Real, power_electric_generated::Real, power_thermal::Real, power_electric_net::Real)
    power_electric_generated = power_electric_generated / 1E6
    power_thermal = power_thermal / 1E6
    power_electric_net = power_electric_net / 1E6
    cost = 27.0 * (land / 1000.0)^0.2 # site
    cost += 111.661 * (building_volume / 80.0e3)^0.62 # tokamak building
    cost += 78.9 * (power_electric_generated / 1246)^0.5 # turbine/generator
    cost += 16.804 * ((power_thermal - power_electric_net) / 1860.0)^0.5 # heat rejection
    cost += 22.878 * (power_electric_net / 1000.0)^0.3 # electrical equipment
    cost += 21.96 * (power_electric_generated / 1246)^0.3  # Plant auxiliary systems
    cost += 4.309 * (power_electric_generated / 1000.0)^0.3 # power core service building
    cost += 1.513 * (power_electric_generated / 1000.0)^0.3  # service water
    cost += 25.0 * (power_thermal / 1759.0)^0.3  # tritium plant and systems
    cost += 7.11 # control room
    cost += 4.7 * (power_electric_net / 1000.0)^0.3  # On-site AC power supply building
    cost += 2.0 # administrative
    cost += 2.0 # site service
    cost += 2.09 # cyrogenic and inert gas storage
    cost += 0.71 # security
    cost += 4.15  # ventilation stack
    return cost
end

"""
    cost_direct_capital_ARIES(::Type{Val{:hot_cell}}, building_volume)

NOTE: https://www.iter.org/mach/HotCell
"""
function cost_direct_capital_ARIES(::Type{Val{:hot_cell}}, building_volume::Real)
    return 0.34 * 111.661 * (building_volume / 80.0e3)^0.62
end

"""
    cost_direct_capital_ARIES(::Type{Val{:heat_transfer_loop_materials}}, power_thermal::Real)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf (warning uses LiPb for blanket)
"""
function cost_direct_capital_ARIES(::Type{Val{:heat_transfer_loop_materials}}, power_thermal::Real) #, costing::IMAS.costing)
    power_thermal = power_thermal / 1E6
    cost = 50.0 * (power_thermal / 2000.0)^0.55 # water
    cost += 125.0 * (power_thermal / 2000.0)^0.55 # LiPb
    cost += 110.0 * (power_thermal / 2000.0)^0.55 # He
    cost += 0.01 * power_thermal # NbIHX
    cost += 50.0 * (power_thermal / 2000.0)^0.55 # Na
    # return rate_translator(cost, 2009, 10thofkind, costing)
    return cost
end

"""
    cost_direct_capital_ARIES(::Type{Val{:balance_of_plant_equipment}}, power_thermal::Real, power_electric_generated::Real)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_direct_capital_ARIES(::Type{Val{:balance_of_plant_equipment}}, power_thermal::Real, power_electric_generated::Real)
    power_thermal = power_thermal / 1E6
    power_electric_generated = power_electric_generated / 1E6
    cost = 350.0 * (power_thermal / 2620.0)^0.7 # Turbine equipment
    cost += 182.98 * (power_electric_generated / 1200.0)^0.5 # Electrical plant equipment
    cost += 87.52 * ((power_thermal - power_electric_generated) / 2300.0) # Heat rejection equipment
    cost += 88.89 * (power_electric_generated / 1200.0)^0.6 # Miscellenous equipment
    return cost
end

"""
    cost_direct_capital_ARIES(::Type{Val{:fuel_cycle_rad_handling}}, power_thermal::Real, power_electric_net::Real)

NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf (warning uses LiPb for blanket)
"""
function cost_direct_capital_ARIES(::Type{Val{:fuel_cycle_rad_handling}}, power_thermal::Real, power_electric_net::Real)
    power_thermal = power_thermal / 1E6
    power_electric_net = power_electric_net / 1E6
    cost = 15.0 * (power_thermal / 1758.0)^0.85 # radioactive material treatment and management
    cost += 70.0 * (power_thermal / 1758.0)^0.8 # Fuel handling and storage
    cost += 100.0 * (power_electric_net / 2000.0)^0.55 # Hot cell maintenance
    cost += 60.0 # Instrumentation and Control
    cost += 8.0 * (power_thermal / 1000.0)^0.8 # Misc power core equipment
    return cost
end

#= ====================== =#
#  yearly operations cost  #
#= ====================== =#
"""
    cost_operations_ARIES(::Type{Val{:operation_maintenance}}, power_electric_generated::Real)

Yearly cost for maintenance [\$M/year]
NOTE: ARIES https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf
"""
function cost_operations_ARIES(::Type{Val{:operation_maintenance}}, power_electric_generated::Real)
    power_electric_generated = power_electric_generated / 1E6
    return 80.0 * (power_electric_generated / 1200.0)^0.5
end

"""
    cost_operations_ARIES(::Type{Val{:tritium_handling}})

Yearly cost for tritium_handling [\$M/year]
!!!!WRONG!!!! Needs estiamte
"""
function cost_operations_ARIES(::Type{Val{:tritium_handling}})
    return 1.0
end

"""
    cost_operations_ARIES(::Type{Val{:blanket_replacement}}, cost_blanket::Real, blanket_lifetime::Real)

Yearly cost for blanket replacement [\$M/year]
"""
function cost_operations_ARIES(::Type{Val{:blanket_replacement}}, cost_blanket::Real, blanket_lifetime::Real)
    return cost_blanket / blanket_lifetime
end

#= =================== =#
#  Decomissioning cost  #
#= =================== =#
"""
    cost_decomissioning_ARIES(::Type{Val{:decom_wild_guess}}, lifetime::Real)

Cost to decommission the plant [\$M]

LIKELY NEEDS FIXING
"""
function cost_decomissioning_ARIES(::Type{Val{:decom_wild_guess}}, lifetime::Real)
    unit_cost = 2.76 # [$M/year]From GASC
    return unit_cost * lifetime
end

function costing_ARIES(dd, par)
    cst = dd.costing
    cost_direct = cst.cost_direct_capital
    cost_ops = cst.cost_operations
    cost_decom = cst.cost_decommissioning

    ###### Direct Capital ######

    ### Tokamak

    # build layers
    sys = resize!(cost_direct.system, "name" => "tokamak")
    for layer in dd.build.layer
        if layer.fs == Int(_lfs_)
            continue # avoid double counting of hfs and lfs layers
        elseif layer.type == Int(_oh_)
            continue # avoid double counting of oh
        end
        c = cost_direct_capital_ARIES(layer)
        if c > 0
            sub = resize!(sys.subsystem, "name" => replace(layer.name, r"^hfs " => ""))
            sub.cost = c
        end
    end

    # PF coils
    for (name, c) in cost_direct_capital_ARIES(dd.pf_active)
        sub = resize!(sys.subsystem, "name" => name)
        sub.cost = c
    end

    # Heating and current drive
    for hcd in vcat(dd.ec_launchers.beam, dd.ic_antennas.antenna, dd.lh_antennas.antenna, dd.nbi.unit)
        c = cost_direct_capital_ARIES(hcd)
        if c > 0
            sub = resize!(sys.subsystem, "name" => uppercase(hcd.name))
            sub.cost = c
        end
    end

    ### Facility

    sys = resize!(cost_direct.system, "name" => "facility")
    
    if ismissing(dd.balance_of_plant.thermal_cycle, :power_electric_generated) || @ddtime(dd.balance_of_plant.power_electric_net) < 0
        @warn("The plant doesn't generate net electricity therefore costing excludes facility estimates")
        power_electric_net = 0.0
        power_thermal = 0.0
        power_electric_generated = 0.0
    else
        power_electric_net = @ddtime(dd.balance_of_plant.power_electric_net) # should be pulse average
        power_thermal = @ddtime(dd.balance_of_plant.thermal_cycle.total_useful_heat_power)
        power_electric_generated = @ddtime(dd.balance_of_plant.thermal_cycle.power_electric_generated)

        for item in vcat(:land, :buildings, :hot_cell, :heat_transfer_loop_materials, :balance_of_plant_equipment, :fuel_cycle_rad_handling)
            sub = resize!(sys.subsystem, "name" => replace(string(item),"_" => " "))
            if item == :land
                sub.cost = cost_direct_capital_ARIES(item, par.land_space, power_electric_generated)
            elseif item == :buildings
                sub.cost = cost_direct_capital_ARIES(item, par.building_volume,
                    par.land_space, power_electric_generated,
                    power_thermal, power_electric_net)
            elseif item == :hot_cell
                sub.cost = cost_direct_capital_ARIES(item, par.building_volume)
            elseif item == :heat_transfer_loop_materials
                sub.cost = cost_direct_capital_ARIES(item, power_thermal)
            elseif item == :balance_of_plant_equipment
                sub.cost = cost_direct_capital_ARIES(item, power_thermal, power_electric_generated)
            elseif item == :fuel_cycle_rad_handling
                sub.cost = cost_direct_capital_ARIES(item, power_thermal, power_electric_net)
            else
                sub.cost = cost_direct_capital_ARIES(item)
            end
        end

    ###### Operations ######
    sys = resize!(cost_ops.system, "name" => "tritium handling")
    sys.yearly_cost = cost_operations_ARIES(:tritium_handling)

    sys = resize!(cost_ops.system, "name" => "maintenance and operators")
    sys.yearly_cost = cost_operations_ARIES(:operation_maintenance, power_electric_generated)

    sys = resize!(cost_ops.system, "name" => "replacements")
    for item in [:blanket_replacement]
        sub = resize!(sys.subsystem, "name" => string(item))
        if item == :blanket_replacement
            tokamak = cost_direct.system[findfirst(system -> system.name == "tokamak", cost_direct.system)]
            blanket_cost = sum(item.cost for item in tokamak.subsystem if item.name == "blanket")
            sub.yearly_cost = cost_operations_ARIES(:blanket_replacement, blanket_cost, par.blanket_lifetime)
        else
            sub.yearly_cost = cost_operations_ARIES(item)
        end
    end
end

    ###### Decomissioning ######
    sys = resize!(cost_decom.system, "name" => "decommissioning")
    sys.cost = cost_decomissioning_ARIES(:decom_wild_guess, par.lifetime)

     ###### Levelized Cost Of Electricity  ###### #ARIES and Sheffield have different formulas for caluclating this so each should live in its respective actor 
    capital_cost_rate = par.interest_rate / (1 - (1 + par.interest_rate)^(-1.0 * par.lifetime))
    lifetime_cost = 0.0
    for year in 1:par.lifetime
        yearly_cost = (capital_cost_rate * cost_direct.cost + cost_ops.yearly_cost + cost_decom.cost / par.lifetime)
        lifetime_cost += (1.0 + par.escalation_fraction) * (1.0 + par.indirect_cost_rate) * yearly_cost
    end
    dd.costing.cost_lifetime = lifetime_cost
    dd.costing.levelized_CoE = (dd.costing.cost_lifetime * 1E6) / (par.lifetime * 24 * 365 * power_electric_net / 1e3 * par.availability)

end
