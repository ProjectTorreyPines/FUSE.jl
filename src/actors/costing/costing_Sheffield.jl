#= ================== =#
#  Dispatch on symbol  #
#= ================== =#
function cost_direct_capital_Sheffield(item::Symbol, args...; kw...)
    return cost_direct_capital_Sheffield(Val{item}, args...; kw...)
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

#= ============== =#
#  materials cost  #
#= ============== =#
#NOTE: material should be priced by Kg
#NOTE: if something is priced by m^3 then it is for a specific part already
function unit_cost(material::AbstractString)
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



#= =================== =#
#  direct capital cost  #
#= =================== =#

#Equation 18 in Generic magnetic fusion reactor revisited, Sheffield and Milora, FS&T 70 (2016)
function cost_direct_capital_Sheffield(::Type{Val{:main_heat_transfer_system}}, power_thermal::Real)
    power_thermal = power_thermal / 1e6 #want this in megawatts 
    cost = 221 * (power_thermal / 4150)^0.6
    return cost
    #return today_dollars(cost)
end

function cost_direct_capital_Sheffield(::Type{Val{:primary_coils}}, bd::IMAS.build)
    primary_coils_hfs = IMAS.get_build(bd, type=IMAS._tf_, fs=_hfs_)
    cost = 1.5 * primary_coils_hfs.volume * unit_cost(primary_coils_hfs.material)
    return cost
end

function cost_direct_capital_Sheffield(::Type{Val{:shielding_gaps}}, bd::IMAS.build)
    shield_hfs = IMAS.get_build(bd, type=IMAS._shield_, fs=_hfs_, return_only_one=false, raise_error_on_missing = false)
    gaps_hfs = IMAS.get_build(bd, type=IMAS._gap_, fs=_hfs_, return_only_one=false, raise_error_on_missing = false)

    cost = 0

    if !ismissing(shield_hfs)
        for layer in shield_hfs
            vol_shield_hfs = layer.volume
            material_shield_hfs = layer.material
            cost += vol_shield_hfs * unit_cost(material_shield_hfs)
        end
    end

    if !ismissing(gaps_hfs)
        for layer in gaps_hfs
            vol_gaps_hfs = layer.volume
            cost += vol_gaps_hfs * 0.29 # $M/m^3 from Table A.VII 
        end 
    end

    return 1.25 * cost
    
    #return today_dollars(1.25 * cost)
end

function cost_direct_capital_Sheffield(::Type{Val{:structure}}, bd::IMAS.build)
    primary_coils_hfs = IMAS.get_build(bd, type=IMAS._tf_, fs=_hfs_)
    cost = 0.75 * primary_coils_hfs.volume * unit_cost("steel")
    return cost
    #return today_dollars(cost)
end

function cost_direct_capital_Sheffield(::Type{Val{:aux_power}}, ec_power::Real, ic_power::Real, lh_power::Real, nbi_power::Real)
    aux_power = ec_power + ic_power + lh_power + nbi_power
    cost = 1.1 * (aux_power * 5.3) * 1e-6 # $5.3 per W in 2017 dollars, 1e-6 to convert to M$ 
    return cost
    #return today_dollars(cost)
end

#Equation 19
function cost_direct_capital_Sheffield(::Type{Val{:balance_of_plant}}, power_electric_net::Real, power_thermal::Real)
    if power_electric_net <= 0
        return 0.0 
    else 
        power_electric_net = power_electric_net / 1e6 #want input for power electric net, power_thermal in megawatts
        power_thermal = power_thermal / 1e6
        cost = 900 + 900 * (power_electric_net / 1200) * (power_thermal / 4150) * 0.6
        #return today_dollars(cost)
        return cost
    end
end

function cost_direct_capital_Sheffield(::Type{Val{:buildings}}, bd::IMAS.build)
    layers = IMAS.get_build(bd, return_only_one = false)
    vol_fusion_island = 0 
    
    for layer in layers
        vol_fusion_island += layer.volume
    end

    cost = 839 * (vol_fusion_island / 5100)^0.67
    #return today_dollars(cost) 
    return cost
end


#= ====================== =#
#      yearly fuel cost    #
#= ====================== =#

#Equation 23 in Generic magnetic fusion reactor revisited, Sheffield and Milora, FS&T 70 (2016) 
function cost_fuel_Sheffield(::Type{Val{:blanket}}, fixed_charge_rate::Real, initial_cost_blanket::Real, availability::Real, lifetime::Real, neutron_flux::Real, blanket_fluence_lifetime::Real, power_electric_net::Real)
#if electric power is generated then add in the blanket replacement cost; if it's not, just return the blanket capital cost
        blanket_capital_cost = 1.1 * initial_cost_blanket * fixed_charge_rate
        blanket_replacement_cost = ((availability * lifetime * neutron_flux / blanket_fluence_lifetime - 1) * initial_cost_blanket) / lifetime #blanket fluence lifetime in MW*yr/m^2
    
        if power_electric_net > 0.0
            return 1.1 * (blanket_capital_cost + blanket_replacement_cost)
        else 
            @warn "Blanket costs do not include replacements since electric power isn't generated"
            return initial_cost_blanket 
        end
    # return today_dollars(cost)
end

#Equation 24
function cost_fuel_Sheffield(::Type{Val{:divertor}}, fixed_charge_rate::Real, initial_cost_divertor::Real, availability::Real, lifetime::Real, thermal_flux::Real, divertor_fluence_lifetime::Real, power_electric_net)
    divertor_capital_cost = 1.1 * initial_cost_divertor * fixed_charge_rate
    divertor_replacement_cost = (availability * lifetime * thermal_flux / divertor_fluence_lifetime - 1) * initial_cost_divertor / lifetime #divertor_lifetime is fluence lifetime so in MW*yr/m^2

    if power_electric_net > 0.0 
        return 1.1 * (divertor_capital_cost + divertor_replacement_cost)
    else 
        @warn "Divertor costs do not include replacements since electric power isn't generated"
        return initial_cost_divertor
    end
    #return today_dollars(cost)
end

#Table III
function cost_fuel_Sheffield(::Type{Val{:aux_power}}, ec_power::Real, ic_power::Real, lh_power::Real, nbi_power::Real)
    return 0.1 * cost_direct_capital_Sheffield(Val{:aux_power}, ec_power, ic_power, lh_power, nbi_power)
    #return today_dollars(cost)
end

function cost_fuel_Sheffield(::Type{Val{:misc_fuel}}, fixed_charge_rate::Real)
    cost = 24 * fixed_charge_rate + 0.4 # $M/yr # 24*FCR for misc. replacements plus 0.4 M$ for fuel costs (in 1983 dollars)
    return cost
    #return today_dollars(cost)
end

#= ====================== =#
#  yearly operations cost  #
#= ====================== =#

function cost_operations_maintenance_Sheffield(power_electric_net)
    power_electric_net = power_electric_net / 1e6 
    return 7.7 * (1200/power_electric_net)^0.5
end

#= =================== =#
#  Decomissioning cost  #
#= =================== =#

function cost_decomissioning_Sheffield(power_electric_net::Real)
    cost = 1.0 # placeholder value, should be 1e-3 dollars per kilowatthour 
    return cost
end

function costing_Sheffield(dd, par)
    cst = dd.costing
    cost_direct = cst.cost_direct_capital
    cost_fuel = cst.cost_fuel
    cost_ops = cst.cost_operations
    cost_decom = cst.cost_decommissioning

    bd = dd.build

    fixed_charge_rate = par.fixed_charge_rate
    availability = par.availability
    lifetime = par.lifetime
    initial_cost_divertor = par.initial_cost_divertor
    initial_cost_blanket = par.initial_cost_blanket
    divertor_fluence_lifetime = par.divertor_fluence_lifetime
    blanket_fluence_lifetime = par.blanket_fluence_lifetime
    construction_lead_time = par.construction_lead_time
    inflate_to_start_year = par.inflate_to_start_year
    construction_start_year = par.construction_start_year
    future_inflation_rate = par.future_inflation_rate

    thermal_flux = 1 # placeholder value for thermal flux on the divertor 

    ec_power = 0
    if !isempty(dd.ec_launchers.beam)
        for num in length(dd.ec_launchers.beam[:])
            ec_power += dd.ec_launchers.beam[num].available_launch_power
        end
    end

    ic_power = 0
    if !isempty(dd.ic_antennas.antenna)
        for num in length(dd.ic_antennas.antenna[:])
            ic_power += dd.ic_antennas.antenna[num].available_launch_power
        end
    end

    lh_power = 0
    if !isempty(dd.lh_antennas.antenna)
        for num in length(dd.lh_antennas.antenna[:])
            lh_power += dd.lh_antennas.antenna[num].available_launch_power
        end
    end

    nbi_power = 0
    if !isempty(dd.nbi.unit)
        for num in length(dd.nbi.unit[:])
            nbi_power += dd.nbi.unit[num].available_launch_power
        end
    end

    flux_r = dd.neutronics.time_slice[].wall_loading.flux_r
    flux_z = dd.neutronics.time_slice[].wall_loading.flux_z
    neutron_flux = sum(sqrt.(flux_r.^2 .+ flux_z.^2) / 1e6)/length(flux_r)  


    if ismissing(dd.balance_of_plant.thermal_cycle, :power_electric_generated) || @ddtime(dd.balance_of_plant.power_electric_net) < 0
        @warn("The plant doesn't generate net electricity therefore costing excludes heat transfer and balance of plant estimates")
        power_electric_net = 0.0
        power_thermal = 0.0
        power_electric_generated = 0.0
    else
        power_electric_net = @ddtime(dd.balance_of_plant.power_electric_net) 
        power_electric_generated = @ddtime(dd.balance_of_plant.thermal_cycle.power_electric_generated)
        power_thermal = @ddtime(dd.balance_of_plant.thermal_cycle.total_useful_heat_power)

    end

    ###### Direct Capital ######
    total_direct_capital_cost = 0
    ##### fusion island
    sys = resize!(cost_direct.system, "name" => "fusion island")

    #main heat transfer system 
    sub = resize!(sys.subsystem, "name" => "main heat transfer system")
    sub.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_direct_capital_Sheffield(:main_heat_transfer_system, power_thermal))
    total_direct_capital_cost += sub.cost

    #primary coils 
    sub = resize!(sys.subsystem, "name" => "primary coils")
    sub.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_direct_capital_Sheffield(:primary_coils, bd))
    total_direct_capital_cost += sub.cost

    #shielding gaps 
    sub = resize!(sys.subsystem, "name" => "shielding gaps")
    sub.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_direct_capital_Sheffield(:shielding_gaps, bd))
    total_direct_capital_cost += sub.cost

    #structure 
    sub = resize!(sys.subsystem, "name" => "structure")
    sub.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_direct_capital_Sheffield(:structure, bd))
    total_direct_capital_cost += sub.cost

    #aux power 
    sub = resize!(sys.subsystem, "name" => "aux power")
    sub.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_direct_capital_Sheffield(:aux_power, ec_power, ic_power, lh_power, nbi_power))
    total_direct_capital_cost += sub.cost

    ##### balance of plant
    sys = resize!(cost_direct.system, "name" => "balance of plant")
    sys.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_direct_capital_Sheffield(:balance_of_plant, power_electric_net, power_thermal))
    total_direct_capital_cost += sys.cost

    ##### buildings 
    sys = resize!(cost_direct.system, "name" => "buildings")
    sys.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_direct_capital_Sheffield(:buildings, bd))
    total_direct_capital_cost += sys.cost

    ###### Fuel ######
    total_fuel_cost = 0

    if power_electric_net > 0.0 #if there is no electric power generated, treat blanket and divertor as direct capital costs instead of fuel costs 
        sys = resize!(cost_fuel.system, "name" => "blanket")
        sys.yearly_cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_fuel_Sheffield(:blanket, fixed_charge_rate, initial_cost_blanket, availability, lifetime, neutron_flux, blanket_fluence_lifetime, power_electric_net))
        total_fuel_cost += sys.yearly_cost
    else 
        sys = resize!(cost_direct.system, "name" => "blanket")
        sys.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_fuel_Sheffield(:blanket, fixed_charge_rate, initial_cost_blanket, availability, lifetime, neutron_flux, blanket_fluence_lifetime, power_electric_net))
        total_direct_capital_cost += sys.cost 
    end

    if power_electric_net > 0.0
        sys = resize!(cost_fuel.system, "name" => "divertor")
        sys.yearly_cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_fuel_Sheffield(:divertor, fixed_charge_rate, initial_cost_divertor, availability, lifetime, thermal_flux, divertor_fluence_lifetime, power_electric_net))
        total_fuel_cost += sys.yearly_cost
    else 
        sys = resize!(cost_direct.system, "name" => "divertor")
        sys.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_fuel_Sheffield(:divertor, fixed_charge_rate, initial_cost_divertor, availability, lifetime, thermal_flux, divertor_fluence_lifetime, power_electric_net))
        total_direct_capital_cost += sys.cost
    end

    sys = resize!(cost_fuel.system, "name" => "aux power")
    sys.yearly_cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_fuel_Sheffield(:aux_power, ec_power, ic_power, lh_power, nbi_power))
    total_fuel_cost += sys.yearly_cost 

    sys = resize!(cost_fuel.system, "name" => "misc. fuel")
    sys.yearly_cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_fuel_Sheffield(:misc_fuel, fixed_charge_rate))
    total_fuel_cost += sys.yearly_cost

    ###### Operations & Maintenance #####
    sys = resize!(cost_ops.system, "name" => "operations and maintenance")
    sys.yearly_cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_operations_maintenance_Sheffield(power_electric_net))
    
    ###### Decomissioning ######
    sys = resize!(cost_decom.system, "name" => "decommissioning")
    sys.cost = future_dollars(inflate_to_start_year, construction_start_year, future_inflation_rate, cost_decomissioning_Sheffield(power_electric_net))

    ###### Levelized Cost Of Electricity  ######
    function indirect_charges(construction_lead_time)
      return 1 + (0.5*construction_lead_time/8)
    end

    function capitalization_factor(construction_lead_time)
      return 1.011^(construction_lead_time + 0.61)
    end     

    total_capital_cost = total_direct_capital_cost * capitalization_factor(construction_lead_time) * indirect_charges(construction_lead_time)

    levelized_CoE = 1e3*(total_capital_cost * fixed_charge_rate + total_fuel_cost + cost_operations_maintenance_Sheffield(power_electric_net)) / (power_electric_net*1e-6 * 8760 * availability) + 0.5e-3
    dd.costing.levelized_CoE = levelized_CoE

end