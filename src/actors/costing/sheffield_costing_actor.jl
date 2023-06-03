#= ===================== =#
#  ActorSheffieldCosting  #
#= ===================== =#

Base.@kwdef mutable struct FUSEparameters__ActorSheffieldCosting{T} <: ParametersActor where {T <: Real}
	_parent::WeakRef = WeakRef(nothing)
	_name::Symbol = :not_set
	construction_lead_time::Entry{T} = Entry{T}("year", "Duration of construction"; default = 8.0)
	fixed_charge_rate::Entry{T} = Entry{T}("-", "Constant dollar fixed charge rate"; default = 0.078)
	initial_cost_blanket::Entry{T} = Entry{T}("\$M", "Cost of initial blanket"; default = 200.0)
	divertor_fluence_lifetime::Entry{T} = Entry{T}("MW*yr/m^2", "Divertor fluence over its lifetime"; default = 10.0)
	blanket_fluence_lifetime::Entry{T} = Entry{T}("MW*yr/m^2", "Blanket fluence over its lifetime"; default = 15.0)
end

mutable struct ActorSheffieldCosting{D,P} <: FacilityAbstractActor
	dd::IMAS.dd{D}
	par::FUSEparameters__ActorSheffieldCosting{P}
	function ActorSheffieldCosting(dd::IMAS.dd{D}, par::FUSEparameters__ActorSheffieldCosting{P}; kw...) where {D<:Real,P<:Real}
		logging_actor_init(ActorSheffieldCosting)
		par = par(kw...)
		return new{D,P}(dd, par)
	end
end

function ActorSheffieldCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
	par = act.ActorSheffieldCosting(kw...)
	actor = ActorSheffieldCosting(dd, par)
	step(actor)
	finalize(actor)
	return actor
end

function _step(actor::ActorSheffieldCosting)
	dd = actor.dd
	par = actor.par
	cst = dd.costing

	cost_direct = cst.cost_direct_capital
	cost_ops = cst.cost_operations

	bd = dd.build

	da = DollarAdjust(dd)

	fixed_charge_rate = par.fixed_charge_rate
	availability = cst.availability
	plant_lifetime = cst.plant_lifetime
	initial_cost_blanket = par.initial_cost_blanket
	divertor_fluence_lifetime = par.divertor_fluence_lifetime
	blanket_fluence_lifetime = par.blanket_fluence_lifetime
	construction_lead_time = par.construction_lead_time

	thermal_flux = 10.0 # placeholder value for thermal flux on the divertor 

	ec_power = 0.0
	if !isempty(dd.ec_launchers.beam)
		for num in length(dd.ec_launchers.beam[:])
			ec_power += dd.ec_launchers.beam[num].available_launch_power
		end
	end

	ic_power = 0.0
	if !isempty(dd.ic_antennas.antenna)
		for num in length(dd.ic_antennas.antenna[:])
			ic_power += dd.ic_antennas.antenna[num].available_launch_power
		end
	end

	lh_power = 0.0
	if !isempty(dd.lh_antennas.antenna)
		for num in length(dd.lh_antennas.antenna[:])
			lh_power += dd.lh_antennas.antenna[num].available_launch_power
		end
	end

	nb_power = 0.0
	if !isempty(dd.nbi.unit)
		for num in length(dd.nbi.unit[:])
			nb_power += dd.nbi.unit[num].available_launch_power
		end
	end

	flux_r = dd.neutronics.time_slice[].wall_loading.flux_r
	flux_z = dd.neutronics.time_slice[].wall_loading.flux_z
	neutron_flux = sum(sqrt.(flux_r .^ 2 .+ flux_z .^ 2) / 1e6) / length(flux_r)

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
	total_direct_capital_cost = 0.0
	##### fusion island
	sys_fi = resize!(cost_direct.system, "name" => "tokamak")

	#main heat transfer system 
	sub = resize!(sys_fi.subsystem, "name" => "main heat transfer system")
	sub.cost = cost_direct_capital_Sheffield(:main_heat_transfer_system, power_thermal, da)
	total_direct_capital_cost += sub.cost

	#primary coils 
	sub = resize!(sys_fi.subsystem, "name" => "primary coils")
	sub.cost = cost_direct_capital_Sheffield(:primary_coils, cst, bd, da)
	total_direct_capital_cost += sub.cost

	#shielding gaps 
	sub = resize!(sys_fi.subsystem, "name" => "shielding gaps")
	sub.cost = cost_direct_capital_Sheffield(:shielding_gaps, cst, bd, da)
	total_direct_capital_cost += sub.cost

	#structure 
	sub = resize!(sys_fi.subsystem, "name" => "structure")
	sub.cost = cost_direct_capital_Sheffield(:structure, cst, bd, da)
	total_direct_capital_cost += sub.cost

	#aux power 
	sub = resize!(sys_fi.subsystem, "name" => "aux power")
	sub.cost = cost_direct_capital_Sheffield(:aux_power, ec_power, ic_power, lh_power, nb_power, da)
	total_direct_capital_cost += sub.cost

	##### Facility structures, buildings and site 
	sys_bld = resize!(cost_direct.system, "name" => "facility")

	#balance of plant
	sub = resize!(sys_bld.subsystem, "name" => "balance of plant")
	sub.cost = cost_direct_capital_Sheffield(:balance_of_plant, power_electric_net, power_thermal, da)
	total_direct_capital_cost += sub.cost

	#buildings 
	sub = resize!(sys_bld.subsystem, "name" => "buildings")
	sub.cost = cost_direct_capital_Sheffield(:buildings, bd, da)
	total_direct_capital_cost += sub.cost

	###### Fuel ######
	total_fuel_cost = 0

	sys = resize!(cost_ops.system, "name" => "blanket")
	sys.yearly_cost = cost_fuel_Sheffield(:blanket, dd, fixed_charge_rate, initial_cost_blanket, availability, plant_lifetime, neutron_flux, blanket_fluence_lifetime, power_electric_net, da)
	total_fuel_cost += sys.yearly_cost

	sys = resize!(cost_ops.system, "name" => "divertor")
	sys.yearly_cost = cost_fuel_Sheffield(:divertor, dd, fixed_charge_rate, availability, plant_lifetime, neutron_flux, thermal_flux, divertor_fluence_lifetime, power_electric_net, da)
	total_fuel_cost += sys.yearly_cost

	sys = resize!(cost_ops.system, "name" => "aux power")
	sys.yearly_cost = cost_fuel_Sheffield(:aux_power, ec_power, ic_power, lh_power, nb_power, da)
	total_fuel_cost += sys.yearly_cost

	sys = resize!(cost_ops.system, "name" => "misc. fuel")
	sys.yearly_cost = cost_fuel_Sheffield(:misc_fuel, fixed_charge_rate, da)
	total_fuel_cost += sys.yearly_cost

	###### Operations & Maintenance #####
	sys = resize!(cost_ops.system, "name" => "operations and maintenance")
	sys.yearly_cost = cost_operations_maintenance_Sheffield(power_electric_net, da)

	###### Levelized Cost Of Electricity  ######
	function indirect_charges(construction_lead_time)
		return 1 + (0.5 * construction_lead_time / 8)
	end

	function capitalization_factor(construction_lead_time)
		return 1.011^(construction_lead_time + 0.61)
	end

	total_capital_cost = total_direct_capital_cost * capitalization_factor(construction_lead_time) * indirect_charges(construction_lead_time)

	cst.levelized_CoE = 1e3 * (total_capital_cost * fixed_charge_rate + total_fuel_cost + cost_operations_maintenance_Sheffield(power_electric_net, da)) / (power_electric_net * 1e-6 * 8760 * availability) + 0.5e-3

	return actor
end

function _finalize(actor::ActorSheffieldCosting)
	# sort system/subsystem by their costs
	sort!(actor.dd.costing.cost_direct_capital.system, by = x -> x.cost, rev = true)
	for sys in actor.dd.costing.cost_direct_capital.system
		sort!(sys.subsystem, by = x -> x.cost, rev = true)
	end

	return actor
end

#= =================== =#
#  direct capital cost  #
#= =================== =#

#Equation 18 in Generic magnetic fusion reactor revisited, Sheffield and Milora, FS&T 70 (2016)
function cost_direct_capital_Sheffield(::Type{Val{:main_heat_transfer_system}}, power_thermal::Real, da::DollarAdjust)
	da.year_assessed = 2016
	power_thermal = power_thermal / 1e6 #want this in megawatts 
	cost = 221 * (power_thermal / 4150)^0.6
	return future_dollars(cost, da)
end

function cost_direct_capital_Sheffield(::Type{Val{:primary_coils}}, cst::IMAS.costing, bd::IMAS.build, da::DollarAdjust)
	da.year_assessed = 2016   # Year the materials costs were assessed 
	primary_coils_hfs = IMAS.get_build_layer(bd.layer, type = IMAS._tf_, fs = _hfs_)
	cost = 1.5 * primary_coils_hfs.volume * unit_cost(bd.tf.technology, cst)
	return future_dollars(cost, da)
end

function cost_direct_capital_Sheffield(::Type{Val{:shielding_gaps}}, cst::IMAS.costing, bd::IMAS.build, da::DollarAdjust)
	da.year_assessed = 2016
	shields_hfs = IMAS.get_build_layers(bd.layer, type = IMAS._shield_, fs = _hfs_)
	gaps_hfs = IMAS.get_build_layers(bd.layer, type = IMAS._gap_, fs = _hfs_)

	cost = 0

    for layer in shields_hfs
        vol_shield_hfs = layer.volume
        material_shield_hfs = layer.material
        cost += vol_shield_hfs * unit_cost(material_shield_hfs, cst)
    end

    for layer in gaps_hfs
        vol_gaps_hfs = layer.volume
        cost += vol_gaps_hfs * 0.29 # $M/m^3 from Table A.VII 
    end

	return future_dollars(1.25 * cost, da)

end

function cost_direct_capital_Sheffield(::Type{Val{:structure}}, cst::IMAS.costing, bd::IMAS.build, da::DollarAdjust)
	da.year_assessed = 2016
	primary_coils_hfs = IMAS.get_build_layer(bd.layer, type = IMAS._tf_, fs = _hfs_)
	cost = 0.75 * primary_coils_hfs.volume * unit_cost("steel", cst)
	return future_dollars(cost, da)
end

function cost_direct_capital_Sheffield(::Type{Val{:aux_power}}, ec_power::Real, ic_power::Real, lh_power::Real, nb_power::Real, da::DollarAdjust)
	da.year_assessed = 2016
	aux_power = ec_power + ic_power + lh_power + nb_power
	cost = 1.1 * (aux_power * 5.3) * 1e-6 # convert to M$ 
	return future_dollars(cost, da)
end

#Equation 19
function cost_direct_capital_Sheffield(::Type{Val{:balance_of_plant}}, power_electric_net::Real, power_thermal::Real, da::DollarAdjust)
	da.year_assessed = 2010  # pg. 19 of Generic magnetic fusion reactor revisited 
	if power_electric_net <= 0
		return 0.0
	else
		power_electric_net = power_electric_net / 1e6 #want input for power electric net, power_thermal in megawatts
		power_thermal = power_thermal / 1e6
		cost = 900 + 900 * (power_electric_net / 1200) * (power_thermal / 4150) * 0.6
		return future_dollars(cost, da)
	end
end

function cost_direct_capital_Sheffield(::Type{Val{:buildings}}, bd::IMAS.build, da::DollarAdjust)
	da.year_assessed = 2016

    vol_fusion_island = 0.0
	for layer in bd.layer
		vol_fusion_island += layer.volume
	end

	cost = 839.0 * (vol_fusion_island / 5100)^0.67
	return future_dollars(cost, da)
end


#= ====================== =#
#      yearly fuel cost    #
#= ====================== =#

#Equation 23 in Generic magnetic fusion reactor revisited, Sheffield and Milora, FS&T 70 (2016) 
function cost_fuel_Sheffield(::Type{Val{:blanket}}, dd::IMAS.dd, fixed_charge_rate::Real, initial_cost_blanket::Real, availability::Real, plant_lifetime::Real, neutron_flux::Real, blanket_fluence_lifetime::Real, power_electric_net::Real, da::DollarAdjust)
	da.year_assessed = 2016

	bd = dd.build
	cst = dd.costing

	blankets = IMAS.get_build_layers(bd.layer, type = IMAS._blanket_, fs = _lfs_)

	if isempty(blankets)
		cost = 0
	else
        blanket = blankets[1]
		initial_cost_blanket = blanket.volume * unit_cost(blanket.material, cst)

		blanket_capital_cost = 1.1 * initial_cost_blanket * fixed_charge_rate
		blanket_replacement_cost = ((availability * plant_lifetime * neutron_flux / blanket_fluence_lifetime - 1) * initial_cost_blanket) / plant_lifetime #blanket fluence lifetime in MW*yr/m^2

		cost = 1.1 * (blanket_capital_cost + blanket_replacement_cost)
	end
	return future_dollars(cost, da)
end

#Equation 24
function cost_fuel_Sheffield(::Type{Val{:divertor}}, dd::IMAS.dd, fixed_charge_rate::Real, availability::Real, plant_lifetime::Real, neutron_flux::Real, thermal_flux::Real, divertor_fluence_lifetime::Real, power_electric_net, da::DollarAdjust)
	da.year_assessed = 2016

	initial_cost_divertor = 0.1 * 0.114 * 0.8 * ((IMAS.fusion_power(dd) / 1e6) / neutron_flux)

	divertor_capital_cost = 1.1 * initial_cost_divertor * fixed_charge_rate
	divertor_replacement_cost = (availability * plant_lifetime * thermal_flux / divertor_fluence_lifetime - 1) * initial_cost_divertor / plant_lifetime #divertor fluence lifetime in MW*yr/m^2

	cost = 1.1 * (divertor_capital_cost + divertor_replacement_cost)
	return future_dollars(cost, da)

end

#Table III
function cost_fuel_Sheffield(::Type{Val{:aux_power}}, ec_power::Real, ic_power::Real, lh_power::Real, nb_power::Real, da::DollarAdjust)
	return 0.1 * cost_direct_capital_Sheffield(Val{:aux_power}, ec_power, ic_power, lh_power, nb_power, da)
end

function cost_fuel_Sheffield(::Type{Val{:misc_fuel}}, fixed_charge_rate::Real, da::DollarAdjust)
	da.year_assessed = 1983
	cost = 24 * fixed_charge_rate + 0.4 # $M/yr # 24*FCR for misc. replacements plus 0.4 M$ for fuel costs (in 1983 dollars)
	return future_dollars(cost, da)
end

#= ====================== =#
#  yearly operations cost  #
#= ====================== =#

function cost_operations_maintenance_Sheffield(power_electric_net::Real, da::DollarAdjust)
	da.year_assessed = 1983
	if power_electric_net == 0
		cost = 49.1 # $M/yr 
	else
		power_electric_net = power_electric_net / 1e6
		cost = 7.7 * (power_electric_net / 1200)^0.5 # there is a typo in Sheffield 1986 equation F.10, which has been corrected here
	end
	return future_dollars(cost, da)
end
