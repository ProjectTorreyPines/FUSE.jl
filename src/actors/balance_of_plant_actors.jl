#= =================== =#
#  ActorBalanceOfPlant  #
#= =================== =#

Base.@kwdef mutable struct ActorBalanceOfPlant <: FacilityAbstractActor
    dd::IMAS.dd
    blanket_multiplier::Real
    efficiency_reclaim::Real
    thermal_electric_conversion_efficiency::Real
end

function ParametersActor(::Type{Val{:ActorBalanceOfPlant}})
    par = ParametersActor(nothing)
    par.blanket_multiplier = Entry(Real, "", "Neutron thermal power multiplier in blanket"; default=1.2)
    par.efficiency_reclaim = Entry(Real, "", "Reclaim efficiency of thermal power hitting the blanket"; default=0.6)
    par.thermal_electric_conversion_efficiency = Entry(Real, "", "Efficiency of the steam cycle, thermal to electric"; default=0.4)
    return par
end

"""
    ActorBalanceOfPlant(dd::IMAS.dd, act::ParametersActor; gasc_method=false, kw...)

Balance of plant actor that estimates the Net electrical power output by estimating the balance of plant electrical needs and compares it to the electricity generated from the thermal cycle.

Setting `gasc_method = true` simply assumes that the power to balance a plant is 7% of the electricity generated.

Setting `gasc_method = false` subdivides the power plant electrical needs to [:cryostat, :tritium_handling, :pumping] using  EU-DEMO numbers.

!!! note 
    Stores data in `dd.balance_of_plant`
"""
function ActorBalanceOfPlant(dd::IMAS.dd, act::ParametersActor; gasc_method=false, kw...)
    par = act.ActorBalanceOfPlant(kw...)
    actor = ActorBalanceOfPlant(dd, par.blanket_multiplier, par.efficiency_reclaim, par.thermal_electric_conversion_efficiency)
    step(actor, gasc_method)
    finalize(actor)
    return actor
end

function step(actor::ActorBalanceOfPlant, gasc_method)
    dd = actor.dd
    bop = dd.balance_of_plant
    empty!(bop)

    bop.time = dd.core_profiles.time
    bop.thermal_cycle.thermal_electric_conversion_efficiency = ones(length(bop.time)) .* actor.thermal_electric_conversion_efficiency

    # ======= #
    # THERMAL #
    # ======= #
    bop_thermal = bop.thermal_cycle

    ### Blanket ###
    sys = resize!(bop_thermal.system, "name" => "blanket", "index" => 1)
    sys.power_in = [sum([bmod.time_slice[time].power_thermal_extracted for bmod in dd.blanket.module]) for time in bop.time]

    ### Divertor ###
    sys = resize!(bop_thermal.system, "name" => "divertors", "index" => 2)
    sys.power_in = sum([IMAS.get_time_array(div.power_thermal_extracted, :data, bop.time, :constant) for div in dd.divertors.divertor])

    # ======== #
    # ELECTRIC #
    # ======== #
    bop_electric = bop.power_electric_plant_operation

    ## H&CD
    sys = resize!(bop_electric.system, "name" => "H&CD", "index" => 1)
    sys.power = zeros(length(bop.time))
    for (idx, hcd_system) in enumerate(intersect([:nbi, :ec_launchers, :ic_antennas, :lh_antennas], keys(dd)))
        sub_sys = resize!(sys.subsystem, "name" => string(hcd_system), "index" => idx)
        sub_sys.power = electricity(getfield(dd, hcd_system), bop.time)
        sys.power .+= sub_sys.power
    end

    ## balance of plant systems
    if gasc_method
        sys = resize!(bop_electric.system, "name" => "BOP_gasc", "index" => 2)
        sys.power = 0.07 .* bop_thermal.power_electric_generated
    else
        # More realistic DEMO numbers
        bop_systems = [:cryostat, :tritium_handling, :pumping, :pf_active] # index 2 : 5
        for (idx, system) in enumerate(bop_systems)
            sys = resize!(bop_electric.system, "name" => string(system), "index" => (idx + 1))
            sys.power = electricity(system, bop.time)
        end
    end
    return actor
end

function heating_and_current_drive_calc(system_unit, time_array::Vector{<:Real})
    power_electric_total = zeros(length(time_array))
    for item_unit in system_unit
        efficiency = prod([getproperty(item_unit.efficiency, i) for i in keys(item_unit.efficiency)])
        power_electric_total .+= IMAS.get_time_array(item_unit.power_launched, :data, time_array, :constant) ./ efficiency
    end
    return power_electric_total
end

function parse_core_sources_sum_heating(cs::IMAS.core_sources, identifier_index::Int64, time_array::Vector{<:Real})
    return IMAS.interp1d(cs.time, [findall(cs.source, "identifier.index" => identifier_index)[1].profiles_1d[t].total_ion_power_inside[end] .+ findall(cs.source, "identifier.index" => 6)[1].profiles_1d[t].electrons.power_inside[end] for t in cs.time], :constant).(time_array)
end

function electricity(symbol::Symbol, time_array::Vector{<:Real})
    return electricity(Val{symbol}, time_array)
end

function electricity(nbi::IMAS.nbi, time_array::Vector{<:Real})
    return heating_and_current_drive_calc(nbi.unit, time_array)
end

function electricity(ec_launchers::IMAS.ec_launchers, time_array::Vector{<:Real})
    return heating_and_current_drive_calc(ec_launchers.launcher, time_array)
end

function electricity(ic_antennas::IMAS.ic_antennas, time_array::Vector{<:Real})
    return heating_and_current_drive_calc(ic_antennas.antenna, time_array)
end

function electricity(lh_antennas::IMAS.lh_antennas, time_array::Vector{<:Real})
    return heating_and_current_drive_calc(lh_antennas.antenna, time_array)
end

# Dummy functions values taken from DEMO 2017  https://iopscience.iop.org/article/10.1088/0029-5515/57/1/016011
function electricity(::Type{Val{:cryostat}}, time_array::Vector{<:Real})
    return 30e6 .* ones(length(time_array)) # MWe
end

function electricity(::Type{Val{:tritium_handling}}, time_array::Vector{<:Real})
    return 15e6 .* ones(length(time_array)) # MWe
end

function electricity(::Type{Val{:pumping}}, time_array::Vector{<:Real})
    return 80e6 .* ones(length(time_array)) # MWe    (Note this should not be a constant!)
end

function electricity(::Type{Val{:pf_active}}, time_array::Vector{<:Real})
    return 0e6 .* ones(length(time_array)) # MWe    (Note this should not be a constant!)
end

function thermal_power(symbol::Symbol, dd::IMAS.dd, actor::ActorBalanceOfPlant, time_array::Vector{<:Real})
    return thermal_power(Val{symbol}, dd, actor, time_array)
end

function thermal_power(::Type{Val{:blanket}}, dd::IMAS.dd, actor::ActorBalanceOfPlant, time_array::Vector{<:Real})
    power_fusion, time_array_fusion = IMAS.total_power_time(dd.core_sources, [6])
    return actor.blanket_multiplier .* IMAS.interp1d(time_array_fusion, 4 .* power_fusion, :constant).(time_array) # blanket_multiplier * P_neutron
end

function thermal_power(::Type{Val{:diverters}}, dd::IMAS.dd, actor::ActorBalanceOfPlant, time_array::Vector{<:Real})
    return actor.efficiency_reclaim .* IMAS.total_power_source(IMAS.total_sources(dd.core_sources, dd.core_profiles.profiles_1d[])) .* ones(length(time_array))
end