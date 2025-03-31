#= ================= =#
#  ActorThermalPlant  #
#= ================= =#
using BalanceOfPlantSurrogate

abstract type AbstractActorThermalPlant{D,P} <: SingleAbstractActor{D,P} end

Base.@kwdef mutable struct FUSEparameters__ActorThermalSystemModels{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

Base.@kwdef mutable struct FUSEparameters__ActorThermalPlant{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:fixed_plant_efficiency, :network, :surrogate], "-", "Power plant heat cycle efficiency"; default=:surrogate)
    fixed_plant_efficiency::Entry{T} = Entry{T}("-", "Overall thermal cycle efficiency (if `model=:fixed_plant_efficiency`)"; default=0.35, check=x -> @assert 1.0 >= x >= 0.0 "must be: 1.0 >= fixed_plant_efficiency >= 0.0")
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorThermalPlant{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorThermalPlant{P}}
    act::ParametersAllActors{P}
    plant_actor::Union{ActorNoOperation{D,P}, AbstractActorThermalPlant{D,P}}
end

function ActorThermalPlant(dd::IMAS.dd{D}, par::FUSEparameters__ActorThermalPlant{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorThermalPlant)
    par = OverrideParameters(par; kw...)
    noop = ActorNoOperation(dd,act.ActorNoOperation)
    return ActorThermalPlant(dd, par, act, noop)
end

"""
    ActorThermalPlant(dd::IMAS.dd, act::ParametersAllActors; kw...)

!!! note

    Provides a common interface to run different thermal plant actors
"""
function ActorThermalPlant(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorThermalPlant(dd, act.ActorThermalPlant, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorThermalPlant)
    dd = actor.dd
    par = actor.par
    act = actor.act

    bop = dd.balance_of_plant

    breeder_heat_load = @ddtime(bop.power_plant.heat_load.breeder)
    divertor_heat_load = @ddtime(bop.power_plant.heat_load.divertor)
    wall_heat_load = @ddtime(bop.power_plant.heat_load.wall)

    # Nothing to do in absence of a blanket
    if breeder_heat_load == 0.0
        empty!(dd.balance_of_plant)
        @warn "No blanket present for ActorThermalPlant to do anything"
        return actor
    end

    # fixed cycle efficiency
    if par.model == :fixed_plant_efficiency
        @ddtime(bop.thermal_efficiency_plant = par.fixed_plant_efficiency)
        @ddtime(bop.power_plant.total_heat_supplied = breeder_heat_load + divertor_heat_load + wall_heat_load)
        @ddtime(bop.power_plant.power_electric_generated = @ddtime(bop.power_plant.total_heat_supplied) * par.fixed_plant_efficiency)
    elseif par.model == :surrogate
        BOP = BalanceOfPlantSurrogate.BOPsurrogate(Symbol(bop.power_plant.power_cycle_type))
        plant_efficiency = BOP(breeder_heat_load, divertor_heat_load, wall_heat_load)
        @ddtime(bop.thermal_efficiency_plant = plant_efficiency)
        @ddtime(bop.power_plant.total_heat_supplied = breeder_heat_load + divertor_heat_load + wall_heat_load)
        @ddtime(bop.power_plant.power_electric_generated = @ddtime(bop.power_plant.total_heat_supplied) * plant_efficiency)
    else
        ext = Base.get_extension(@__MODULE__, :ThermalSystemModelsExt)
        if ext === nothing
            throw(MissingExtensionError("ActorThermalSystemModels", "ThermalSystemModels"))
        end
        actor.plant_actor = ext.ActorThermalSystemModels(dd, act.ActorThermalSystemModels; par.verbose, par.do_plot)
        step(actor.plant_actor)
    end

    return actor
end

function _finalize(actor::ActorThermalPlant)
    par = actor.par

    if par.model == :fixed_plant_efficiency
        #pass
    elseif par.model == :surrogate
        #pass
    else
        finalize(actor.plant_actor)
    end

    return actor
end