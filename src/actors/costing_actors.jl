#= ============== =#
#  materials cost #
#= ============== =#
#NOTE: material should be priced by Kg
#NOTE: if something is priced by m^3 then it is for a specific part already
function unit_cost(material::String)
    if material == "Vacuum"
        return 0.0 # $M/m^3
    elseif material == "ReBCO"
        return 87.5 / 2 # $M/m^3
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

function cost_direct_capital(layer::IMAS.build__layer)
    if layer.type == Int(_oh_)
        build = IMAS.parent(IMAS.parent(layer))
        return unit_cost(build.oh.technology) * layer.volume
    elseif layer.type == Int(_tf_)
        build = IMAS.parent(IMAS.parent(layer))
        return unit_cost(build.tf.technology) * layer.volume
    elseif layer.type == Int(_shield_)
        return layer.volume * 0.29  # $M/m^3
    elseif layer.type == Int(_blanket_)
        return layer.volume * 0.75  # $M/m^3
    elseif layer.type ∈ [Int(_wall_), Int(_vessel_), Int(_cryostat_)]
        return layer.volume * 0.36  # $M/m^3
    else
        return unit_cost(layer.material) * layer.volume
    end
end

function cost_direct_capital(ecl::IMAS.ec_launchers__beam)
    ecl.available_launch_power / 1E6 * 3.0 # $/W #ARIES
end

function cost_direct_capital(ica::IMAS.ic_antennas__antenna)
    ica.available_launch_power / 1E6 * 1.64 #$/W ARIES
end

function cost_direct_capital(lha::IMAS.lh_antennas__antenna)
    lha.available_launch_power / 1E6 * 2.13 #$/W ARIES
end

function cost_direct_capital(nbu::IMAS.nbi__unit)
    nbu.available_launch_power / 1E6 * 4.93 #$/W ARIES
end

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

function cost_direct_capital(pf_active::IMAS.pf_active)
    dd = IMAS.top_dd(pf_active)
    c = Dict("OH" => 0.0, "PF" => 0.0)
    for coil in pf_active.coil
        if coil.name == "OH"
            c["OH"] += cost_direct_capital(coil, dd.build.oh.technology)
        else
            c["PF"] += cost_direct_capital(coil, dd.build.pf_active.technology)
        end
    end
    return c
end

function cost_direct_capital(coil::IMAS.pf_active__coil, technology::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology})
    return IMAS.volume(coil) * unit_cost(technology)
end

# NOTE costs are based on ARES study https://cer.ucsd.edu/_files/publications/UCSD-CER-13-01.pdf

function cost_direct_capital(::Type{Val{:land}}, land::Real, power_electric_net::Real)
    1.2 * land * 20.0e-3 * (power_electric_net / 1000.0)^0.3
end
cost_direct_capital(item::Symbol, land, power_electric_generated) = cost_direct_capital(Val{item}, land, power_electric_generated)

function cost_direct_capital(::Type{Val{:buildings}}, land::Real, building_volume::Real, power_electric_generated::Real, power_thermal::Real, power_electric_net::Real) # ARIES
    cost = 27.0 * (land / 1000.0)^0.2 # site
    cost += 111.661 * (building_volume / 80.0e3)^0.62 # tokamak building
    cost += 78.9 * (power_electric_generated / 1246)^0.5 # turbine/generator
    cost += 16.804 * ((power_thermal - power_electric_net) / 1860.0)^0.5 # heat rejection
    cost += 22.878 * (power_electric_net / 1000.0)^0.3 # electrical equipment
    cost += 21.96 * (power_electric_generated / 1246)^0.3  # Plant Auxiliary Systems
    cost += 4.309 * (power_electric_generated / 1000.0)^0.3 # power core service building
    cost += 1.513 * (power_electric_generated / 1000.0)^0.3  # service water
    cost += 25.0 * (power_thermal / 1759.0)^0.3  # fuel handling
    cost += 7.11 # control room
    cost += 4.7 * (power_electric_net / 1000.0)^0.3  # On-site AC Power Supply Building
    cost += 2.0 # administrative
    cost += 2.0 # site service
    cost += 2.09 # cyrogenic and inert gas storage
    cost += 0.71 # security
    cost += 4.15  # Ventilation Stack
    return cost
end

cost_direct_capital(item::Symbol,
    building_volume,
    land,
    power_electric_generated,
    power_thermal,
    power_electric_net) = cost_direct_capital(
    Val{item},
    building_volume,
    land,
    power_electric_generated,
    power_thermal,
    power_electric_net)

cost_direct_capital(item::Symbol, power_electric_generated) = cost_direct_capital(Val{item}, power_electric_generated)

function cost_direct_capital(::Type{Val{:hot_cell}}, building_volume) # https://www.iter.org/mach/HotCell
    0.34 * 111.661 * (building_volume / 80.0e3)^0.62
end

function cost_direct_capital(::Type{Val{:heat_transfer_loop_materials}},power_thermal) # ARIES (warning uses LiPb for blanket)
    cost = 50 * (power_thermal / 2000.0) ^ 0.55 # water
    cost += 125 * (power_thermal / 2000.0) ^ 0.55 # LiPb
    cost +=  110.0 * (power_thermal / 2000.0) ^ 0.55 # He
    cost += 0.01 * power_thermal # NbIHX
    cost += 50.0 * (power_thermal / 2000.0) ^ 0.55 # Na
    return cost
end

function cost_direct_capital(::Type{Val{:balance_of_plant_equipment}},power_thermal, power_electric_generated) # ARIES
    cost = 350 * (power_thermal / 2620.0) ^ 0.7 # Turbine equipment
    cost += 182.98 * (power_electric_generated / 1200.0) ^ 0.5 # Electrical plant equipment
    cost +=  87.52 * ((power_thermal - power_electric_generated) / 2300.0) # Heat rejection equipment
    cost += 88.89 * (power_electric_generated / 1200.0) ^ 0.6 # Miscellenous equipment
    return cost
end

function cost_direct_capital(::Type{Val{:fuel_cycle_rad_handling}},power_thermal, power_electric_net) # ARIES (warning uses LiPb for blanket)
    cost = 15 * (power_thermal / 1758.0) ^ 0.85 # radioactive material treatment and management
    cost += 70 * (power_thermal / 1758.0) ^ 0.8 # Fuel handling and storage
    cost +=  100. * (power_electric_net / 2000.0) ^ 0.55 # Hot cell maintanance
    cost += 60. # Instrumentation and Control
    cost += 8 * (power_thermal / 1000.0) ^ 0.8 # Misc power core equipment
    return cost
end

function cost_operations(::Type{Val{:operation_maintanance}}, power_electric_generated)
    80.0 * (power_electric_generated / 1200.0)^0.5
end

function cost_operations(::Type{Val{:fuel}})
    1.0
end
cost_operations(item::Symbol) = cost_operations(Val{item})
cost_operations(item::Symbol, a) = cost_operations(Val{item}, a)
cost_operations(item::Symbol, power_electric_generated, cost_blanket) = cost_operations(Val{item}, power_electric_generated, cost_blanket)

function cost_operations(::Type{Val{:blanket_replacement}}, cost_blanket, blanket_lifetime) # find blanket and replace every x-years
    cost_blanket * 1.2 / blanket_lifetime # 20% contingency for extras
end

function cost_decomissioning(::Type{Val{:decom_wild_guess}})
    2.76 # gasc comment needs revisiting
end
cost_decomissioning(item::Symbol) = cost_decomissioning(Val{item})

#= ============ =#
#  ActorCosting  #
#= ============ =#

mutable struct ActorCosting <: FacilityAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    function ActorCosting(dd::IMAS.dd, par::ParametersActor; kw...)
        par = par(kw...)
        return new(dd, par)
    end
end

function ParametersActor(::Type{Val{:ActorCosting}})
    par = ParametersActor(nothing)
    par.land_space = Entry(Real, "acres", "Plant site space required in acres"; default=1000.0)
    par.building_volume = Entry(Real, "m^3", "Volume of the tokmak building"; default=140.0e3)
    par.intrest_rate = Entry(Real, "", "Anual intrest rate fraction of direct capital cost"; default=0.05)
    par.indirect_cost_rate = Entry(Real, "", "Indirect cost associated with construction, equipment, services, energineering construction management and owners cost"; default=0.4)
    par.lifetime = Entry(Real, "years", "lifetime of the plant"; default=40)
    par.availability = Entry(Real, "", "availability fraction of the plant"; default=0.803)
    par.escalation_fraction = Entry(Real, "", "yearly escalation fraction based on risk assessment"; default=0.05)
    par.blanket_lifetime = Entry(Real, "years", "lifetime of the blanket"; default=6.8)
    return par
end

"""
    ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor estimates the cost of the fusion power plant.

!!! note 
    Stores data in `dd.costing`
"""
function ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCosting(kw...)
    actor = ActorCosting(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function step(actor::ActorCosting)
    par = actor.par
    dd = actor.dd
    cst = dd.costing
    cost_direct = cst.cost_direct_capital
    cost_ops = cst.cost_operations
    cost_decom = cst.cost_decommissioning

    ###### Direct Capital ######

    empty!(cost_direct)

    ### Tokamak

    # build layers
    sys = resize!(cost_direct.system, "name" => "tokamak")
    for layer in dd.build.layer
        if layer.fs == Int(_lfs_)
            continue # avoid double counting of hfs and lfs layers
        elseif layer.type == Int(_oh_)
            continue # avoid double counting of oh
        end
        c = cost_direct_capital(layer)
        if c > 0
            sub = resize!(sys.subsystem, "name" => replace(layer.name, r"^hfs " => ""))
            sub.cost = c
        end
    end

    # PF coils
    for (name, c) in cost_direct_capital(dd.pf_active)
        sub = resize!(sys.subsystem, "name" => name)
        sub.cost = c
    end

    # Heating and current drive
    for hcd in vcat(dd.ec_launchers.beam, dd.ic_antennas.antenna, dd.lh_antennas.antenna, dd.nbi.unit)
        c = cost_direct_capital(hcd)
        if c > 0
            sub = resize!(sys.subsystem, "name" => hcd.name)
            sub.cost = c
        end
    end

    ### Facility
    sys = resize!(cost_direct.system, "name" => "Facility structures, buildings and site")

    if @ddtime(dd.balance_of_plant.thermal_cycle.power_electric_generated) < 0
        @warn("The plant doesn't generate net electricity therefore costing excludes facility estimates")
    else
        power_electric_net = @ddtime(dd.balance_of_plant.power_electric_net) / 1e6 # should be pulse average
        display(dd.balance_of_plant.thermal_cycle)
        power_thermal = @ddtime(dd.balance_of_plant.thermal_cycle.power_thermal_convertable_total) / 1e6 
        power_electric_generated = @ddtime(dd.balance_of_plant.thermal_cycle.power_electric_generated) / 1e6
        for item in vcat(:land, :buildings, :hot_cell, :heat_transfer_loop_materials, :balance_of_plant_equipment, :fuel_cycle_rad_handling)
            sub = resize!(sys.subsystem, "name" => string(item))
            if item == :land
                sub.cost = cost_direct_capital(item, par.land_space, power_electric_generated)
            elseif item == :buildings
                sub.cost = cost_direct_capital(item, par.building_volume,
                    par.land_space, power_electric_generated,
                    power_thermal, power_electric_net)
            elseif item == :hot_cell
                sub.cost = cost_direct_capital(item, par.building_volume)
            elseif item == :heat_transfer_loop_materials
                sub.cost = cost_direct_capital(item, power_thermal)
            elseif item == :balance_of_plant_equipment
                sub.cost = cost_direct_capital(item, power_thermal, power_electric_generated)
            elseif item == :fuel_cycle_rad_handling
                sub.cost = cost_direct_capital(item, power_thermal, power_electric_net)
            else
                sub.cost = cost_direct_capital(item)
            end
        end
    end

    ###### Operations ######
    empty!(cost_ops)

    blanket_cost = sum([item.cost for item in cost_direct.system[1].subsystem if item.name == "blanket"]) # system[1] is always tokamak
    blanket_cost = 275.2622
    sys = resize!(cost_ops.system, "name" => "fuel cycle")
    sys.cost = cost_operations(:fuel)

    sys = resize!(cost_ops.system, "name" => "maintanance and operators")
    sys.cost = cost_operations(:operation_maintanance, power_electric_generated)

    sys = resize!(cost_ops.system, "name" => "replacements")
    for item in [:blanket_replacement]
        sub = resize!(sys.subsystem, "name" => string(item))
        if item == :blanket_replacement
            sub.cost = cost_operations(:blanket_replacement, blanket_cost, par.blanket_lifetime)
        else
            sub.cost = cost_operations(item)
        end
    end

    ###### Decomissioning ######
    empty!(cost_decom)

    sys = resize!(cost_decom.system, "name" => "decommissioning guess")
    sys.cost = cost_decomissioning(:decom_wild_guess)

    ### Levelized Cost Of Electricity 
    capital_cost_rate = par.intrest_rate / (1 - (1 + par.intrest_rate)^(-1.0 * par.lifetime))
    total_cost = 0.0
    for year in 1:par.lifetime
        yearly_cost = (capital_cost_rate * cost_direct.cost + cost_ops.cost + cost_decom.cost)
        total_cost += (1. + par.escalation_fraction) * (1. + par.indirect_cost_rate) * yearly_cost
    end
    dd.costing.levelized_CoE = total_cost / (par.lifetime * 8760 * power_electric_net / 1e3 * par.availability)
    dd.costing.cost_lifetime = total_cost
    return actor
end

function finalize(actor::ActorCosting)
    # sort system/subsystem costs
    sort!(actor.dd.costing.cost_direct_capital.system, by=x -> x.cost, rev=true)
    for sys in actor.dd.costing.cost_direct_capital.system
        sort!(sys.subsystem, by=x -> x.cost, rev=true)
    end
end