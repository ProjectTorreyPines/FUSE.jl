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

function cost_direct_capital(::Type{Val{:land}}, land::Real, power_electric_net::Real)
    1.2 * land * 20.0e-3 * (power_electric_net / 1000.0) ^ 0.3 
end
cost_direct_capital(item::Symbol, land, power_electric_net) = cost_direct_capital(Val{item}, land, power_electric_net)

function cost_direct_capital(::Type{Val{:buildings}},land::Real, building_volume::Real, power_electric_net::Real, power_thermal::Real) # ARIES
    cost = 27.0 * (land / 1000.0)^0.2
    cost += 111.661 * (building_volume / 80.0e3) ^ 0.62 # tokamak building
    cost += 4.309 * (power_electric_net / 1000.0) ^ 0.3 # power core service building
    cost +=  1.513 * (power_electric_net / 1000.0) ^ 0.3  # service water
    cost +=  25.0 * (power_thermal / 1759.0) ^ 0.3  # fuel handling
    cost += 7.11 # control room
    cost += 2.0 # site service
    cost += 2.0 # administrative
    cost += 2.09 # cyrogenic and inert gas storage
    cost += 0.71 # security
    cost += 22.878 * (power_electric_net / 1000.0) ^ 0.3 # service building
    cost += 4.7 * (power_electric_net / 1000.0) ^ 0.3 + 4.15 # On-site AC Power Supply and ventilation    
    return cost
end
cost_direct_capital(item::Symbol,  building_volume,land, power_electric_net, power_thermal) = cost_direct_capital(Val{item},  building_volume, land, power_electric_net, power_thermal)

function cost_direct_capital(::Type{Val{:turbine}},power_electric_generated::Real)
    78.9 * (power_electric_generated / 1246) ^ 0.5 
end
cost_direct_capital(item::Symbol, power_electric_generated) = cost_direct_capital(Val{item}, power_electric_generated)

function cost_direct_capital(::Type{Val{:heat_rejection}}, power_electric_net, power_thermal)
    16.804 * ((power_thermal - power_electric_net) / 1860.0) ^ 0.5
end
cost_direct_capital(item::Symbol, power_electric_net, power_thermal) = cost_direct_capital(Val{item}, power_electric_net, power_thermal)

function cost_direct_capital(::Type{Val{:electrical_equipment}}, power_electric_net)
    22.878 * (power_electric_net / 1000.0) ^ 0.3
end
cost_direct_capital(item::Symbol, power_electric_net) = cost_direct_capital(Val{item}, power_electric_net)

function cost_operations(::Type{Val{:operation_maintanance}}, power_electric_net)
    80.0 * (power_electric_net / 1200.0) ^ 0.5
end
cost_operations(item::Symbol, power_electric_net) = cost_operations(Val{item}, power_electric_net)

function cost_operations(::Type{Val{:fuel}})
    1.0
end
cost_operations(item::Symbol, power_electric_net) = cost_operations(Val{item}, power_electric_net)

function cost_operations(::Type{Val{:blanket_replacement}},cost_blanket) # find blanket and replace every x-years
    cost_blanket * 1.2
end
cost_operations(item::Symbol, cost_blanket) = cost_operations(Val{item}, cost_blanket)

function cost_decomissioning(::Type{Val{:hot_cell}},building_volume) # https://www.iter.org/mach/HotCell
    0.4 * 111.661 * (building_volume / 80.0e3) ^ 0.62
end

function cost_decomissioning(::Type{Val{:decom_wild_guess}})
    2.76e6 # gasc comment needs revisiting
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
    par.land_space = Entry(Real, "m^2","Plant site space required in m²";default=1000.) 
    par.building_volume = Entry(Real, "m^3", "Volume of the tokmak building"; default=140.0e3)
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
    sys = resize!(cost_direct.system, "name" => "facility")

    if @ddtime(dd.balance_of_plant.power_electric_net) < 0
        @warn("The plant doesn't generate net electricity therefore costing excludes facility estimates")
    else
        power_electric_net = @ddtime(dd.balance_of_plant.power_electric_net) / 1e6 # should be pulse average
        power_thermal = sum([maximum(sys.power_in) for sys in dd.balance_of_plant.thermal_cycle.system]) / 1e6 # should be pulse average
        power_electric_generated = @ddtime(dd.balance_of_plant.thermal_cycle.power_electric_generated) / 1e6
        for item in vcat(:land, :buildings, :turbine, :heat_rejection, :electrical_equipment)
            sub = resize!(sys.subsystem, "name" => string(item))
            if item == :land
                sub.cost = cost_direct_capital(item, par.land_space, power_electric_net)
            elseif item == :buildings
                sub.cost = cost_direct_capital(item, par.building_volume, par.land_space, power_electric_net, power_thermal)
            elseif item == :turbine
                sub.cost = cost_direct_capital(item, power_electric_generated)
            elseif item == :heat_rejection
                sub.cost = cost_direct_capital(item, power_electric_net, power_thermal)
            elseif item == :electrical_equipment
                sub.cost = cost_direct_capital(item, power_electric_net)
            else
                sub.cost = cost_direct_capital(item)
            end
        end
    end

    """
    # Fuel Cycle
    sys = resize!(cost_ops.system, "name" => "fuel cycle")

    ###### Decomissioning ######

    # Radioactive waste treatment?
    sys = resize!(cost_decom.system, "name" => "radioactive waste treatment")

    # Demolition
    sys = resize!(cost_decom.system, "name" => "demolition")
    """

    return actor
end

function finalize(actor::ActorCosting)
    # sort system/subsystem costs
    sort!(actor.dd.costing.cost_direct_capital.system, by=x -> x.cost, rev=true)
    for sys in actor.dd.costing.cost_direct_capital.system
        sort!(sys.subsystem, by=x -> x.cost, rev=true)
    end
end