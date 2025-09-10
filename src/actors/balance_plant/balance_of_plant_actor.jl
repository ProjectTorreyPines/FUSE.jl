#= =================== =#
#  ActorBalanceOfPlant  #
#= =================== =#
@actor_parameters_struct ActorBalanceOfPlant{T} begin
end

mutable struct ActorBalanceOfPlant{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorBalanceOfPlant{P}}
    act::ParametersAllActors{P}
    thermal_plant_actor::ActorThermalPlant{D}
    power_needs_actor::ActorPowerNeeds{D}
end

"""
    ActorBalanceOfPlant(dd::IMAS.dd, act::ParametersAllActors; kw...)

Balance of plant actor that estimates the net electrical power output by comparing the balance of plant electrical needs with the electricity generated from the thermal cycle.

!!! note

    Stores data in `dd.balance_of_plant`
"""
function ActorBalanceOfPlant(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorBalanceOfPlant(dd, act.ActorBalanceOfPlant, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorBalanceOfPlant(dd::IMAS.dd, par::FUSEparameters__ActorBalanceOfPlant, act::ParametersAllActors; kw...)
    logging_actor_init(ActorBalanceOfPlant)
    par = OverrideParameters(par; kw...)

    # set the time
    @ddtime(dd.balance_of_plant.time = dd.global_time)

    # set the heat sources
    bop = dd.balance_of_plant
    breeder_heat_load = isempty(dd.blanket.module) ? 0.0 : sum(bmod.time_slice[].power_thermal_extracted for bmod in dd.blanket.module)
    @ddtime(bop.power_plant.heat_load.breeder = breeder_heat_load)
    divertor_heat_load = isempty(dd.divertors.divertor) ? 0.0 : sum((@ddtime(div.power_incident.data)) for div in dd.divertors.divertor)
    @ddtime(bop.power_plant.heat_load.divertor = divertor_heat_load)
    wall_heat_load = abs(IMAS.radiation_losses(dd.core_sources))
    @ddtime(bop.power_plant.heat_load.wall = wall_heat_load)

    # setup actors
    thermal_plant_actor = ActorThermalPlant(dd, act.ActorThermalPlant, act)
    power_needs_actor = ActorPowerNeeds(dd, act.ActorPowerNeeds)
    return ActorBalanceOfPlant(dd, par, act, thermal_plant_actor, power_needs_actor)
end

function _step(actor::ActorBalanceOfPlant)
    finalize(step(actor.thermal_plant_actor))
    finalize(step(actor.power_needs_actor))
    return actor
end

