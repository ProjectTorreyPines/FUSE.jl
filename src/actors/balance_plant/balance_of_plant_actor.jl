#= =================== =#
#  ActorBalanceOfPlant  #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorBalanceOfPlant{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol   = :not_set
    generator_conversion_efficiency::Entry{T} = Entry{T}("-", "Efficiency of the generator"; default=0.95) #  Appl. Therm. Eng. 76 (2015) 123–133, https://doi.org/10.1016/j.applthermaleng.2014.10.093
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
end

mutable struct ActorBalanceOfPlant{D,P} <: FacilityAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorBalanceOfPlant{P}
    act::ParametersAllActors
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
    par = par(kw...)

    # set the time
    @ddtime(dd.balance_of_plant.time = dd.global_time)

    thermal_plant_actor = ActorThermalPlant(dd, act.ActorThermalPlant)
    power_needs_actor   = ActorPowerNeeds(dd, act.ActorPowerNeeds)
    return ActorBalanceOfPlant(dd, par, act, thermal_plant_actor, power_needs_actor)
end

function _step(actor::ActorBalanceOfPlant; kw...)
    dd  = actor.dd
    par = actor.par
    bop = dd.balance_of_plant

    bop_thermal = bop.power_plant
    @ddtime(bop_thermal.generator_conversion_efficiency = par.generator_conversion_efficiency)

    finalize(step(actor.thermal_plant_actor; doplot = par.do_plot))
    finalize(step(actor.power_needs_actor))

    return actor
end
