#= ===================== =#
#  ActorCostingFUSE  #
#= ===================== =#

Base.@kwdef mutable struct FUSEparameters__ActorCostingFUSE{T<:Real} <: ParametersActorBuild{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    construction_lead_time::Entry{T} = Entry{T}("year", "Duration of construction"; default=8.0)
    fixed_charge_rate::Entry{T} = Entry{T}("-", "Constant dollar fixed charge rate"; default=0.078)
    capitalize_blanket::Entry{Bool} = Entry{Bool}("-", "If true, include cost of 1st blanket in direct captial cost"; default=true)
    capitalize_divertor::Entry{Bool} = Entry{Bool}("-", "If true, include cost of 1st divertor in direct captial cost"; default=true)
    divertor_fluence_lifetime::Entry{T} = Entry{T}("MW*yr/m²", "Divertor fluence over its lifetime"; default=10.0)
    blanket_fluence_lifetime::Entry{T} = Entry{T}("MW*yr/m²", "Blanket fluence over its lifetime"; default=15.0)
end

mutable struct ActorCostingFUSE{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCostingFUSE{P}
    function ActorCostingFUSE(dd::IMAS.dd{D}, par::FUSEparameters__ActorCostingFUSE{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorCostingFUSE)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorCostingFUSE(dd::IMAS.dd, act::ParametersAllActors; kw...)

!!! note

    Stores data in `dd.costing`
"""
function ActorCostingFUSE(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCostingFUSE(kw...)
    actor = ActorCostingFUSE(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorCostingFUSE)
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
    divertor_fluence_lifetime = par.divertor_fluence_lifetime
    blanket_fluence_lifetime = par.blanket_fluence_lifetime
    construction_lead_time = par.construction_lead_time

    thermal_flux = 10.0 # placeholder value for thermal flux on the divertor 

    ec_power = 0.0
    if !isempty(dd.ec_launchers.beam)
        for num in eachindex(dd.ec_launchers.beam)
            ec_power += dd.ec_launchers.beam[num].available_launch_power
        end
    end

    ic_power = 0.0
    if !isempty(dd.ic_antennas.antenna)
        for num in eachindex(dd.ic_antennas.antenna)
            ic_power += dd.ic_antennas.antenna[num].available_launch_power
        end
    end

    lh_power = 0.0
    if !isempty(dd.lh_antennas.antenna)
        for num in eachindex(dd.lh_antennas.antenna)
            lh_power += dd.lh_antennas.antenna[num].available_launch_power
        end
    end

    nb_power = 0.0
    if !isempty(dd.nbi.unit)
        for num in eachindex(dd.nbi.unit)
            nb_power += dd.nbi.unit[num].available_launch_power
        end
    end

    wall_loading = dd.neutronics.time_slice[].wall_loading
    flux_r = wall_loading.flux_r
    flux_z = wall_loading.flux_z
    neutron_flux = sum(sqrt.(flux_r .^ 2 .+ flux_z .^ 2) / 1e6) / length(flux_r)

    power_thermal, power_electric_generated, power_electric_net = bop_powers(dd.balance_of_plant)

    ###### Direct Capital ######
    total_direct_capital_cost = 0.0
    ##### fusion island
    sys_fi = resize!(cost_direct.system, "name" => "tokamak")

    # main heat transfer system 
    sub = resize!(sys_fi.subsystem, "name" => "main heat transfer system")
    sub.cost = cost_direct_capital_Sheffield(:main_heat_transfer_system, power_thermal, da)
    total_direct_capital_cost += sub.cost

    # primary coils 
    sub = resize!(sys_fi.subsystem, "name" => "primary coils")
    sub.cost = cost_direct_capital_Sheffield(:primary_coils, cst, bd, da)
    total_direct_capital_cost += sub.cost

    # shields 
    sub = resize!(sys_fi.subsystem, "name" => "shields")
    sub.cost = cost_direct_capital_Sheffield(:shields, cst, bd, da)
    total_direct_capital_cost += sub.cost

    # structure 
    sub = resize!(sys_fi.subsystem, "name" => "structure")
    sub.cost = cost_direct_capital_Sheffield(:structure, cst, bd, da)
    total_direct_capital_cost += sub.cost

    # aux power 
    sub = resize!(sys_fi.subsystem, "name" => "aux power")
    sub.cost = cost_direct_capital_Sheffield(:aux_power, ec_power, ic_power, lh_power, nb_power, da)
    total_direct_capital_cost += sub.cost

    # blanket
    sub = resize!(sys_fi.subsystem, "name" => "blanket")
    sub.cost = cost_direct_capital_Sheffield(:blanket, par.capitalize_blanket, dd, da)
    total_direct_capital_cost += sub.cost

    # divertor
    sub = resize!(sys_fi.subsystem, "name" => "divertor")
    sub.cost = cost_direct_capital_Sheffield(:divertor, par.capitalize_divertor, dd, da)
    total_direct_capital_cost += sub.cost

    ##### Facility structures, buildings and site 
    sys_bld = resize!(cost_direct.system, "name" => "facility")

    # balance of plant
    sub = resize!(sys_bld.subsystem, "name" => "balance of plant")
    sub.cost = cost_direct_capital_FUSE(:balance_of_plant, # whatever your new function arguments are + da)
    total_direct_capital_cost += sub.cost

    # buildings 
    sub = resize!(sys_bld.subsystem, "name" => "buildings")
    sub.cost = cost_direct_capital_Sheffield(:buildings, bd, da)
    total_direct_capital_cost += sub.cost

    ###### Fuel ######
    total_fuel_cost = 0.0

    sys = resize!(cost_ops.system, "name" => "blanket")
    sys.yearly_cost =
        cost_fuel_Sheffield(:blanket, par.capitalize_blanket, dd, fixed_charge_rate, availability, plant_lifetime, neutron_flux, blanket_fluence_lifetime, power_electric_net, da)
    total_fuel_cost += sys.yearly_cost

    sys = resize!(cost_ops.system, "name" => "divertor")
    sys.yearly_cost = cost_fuel_Sheffield(
        :divertor,
        par.capitalize_divertor,
        dd,
        fixed_charge_rate,
        availability,
        plant_lifetime,
        neutron_flux,
        thermal_flux,
        divertor_fluence_lifetime,
        power_electric_net,
        da
    )
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

    cst.levelized_CoE =
        1e3 * (total_capital_cost * fixed_charge_rate + total_fuel_cost + cost_operations_maintenance_Sheffield(power_electric_net, da)) /
        (power_electric_net * 1e-6 * 8760 * availability) + 0.5e-3

    return actor
end

function _finalize(actor::ActorCostingFUSE)
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

function cost_direct_capital_FUSE(::Type{Val{:balance_of_plant}}, # function parameters here, making sure to include da::DollarAdjust)
    da.year_assessed = # put the year that the cost was assessed 

    cost = # some scaling law or laws that define the cost of the components 
    return future_dollars(cost, da)
end