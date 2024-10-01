#= ================= =#
#  ActorThermalPlant  #
#= ================= =#

Base.@kwdef mutable struct FUSEparameters__ActorThermalPlant{T<:Real} <: ParametersActorBuild{T}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:fixed_cycle_efficiency, :network], "-", "Power plant heat cycle efficiency"; default=:network)
    fixed_cycle_efficiency::Entry{T} = Entry{T}("-", "Overall thermal cycle efficiency (if `model=:fixed_cycle_efficiency`)"; default=0.35, check=x -> @assert 1.0 >= x >= 0.0 "must be: 1.0 >= rho_0 >= 0.0")
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

abstract type AbstractActorThermalPlant{D,P} <: SingleAbstractActor{D,P} end

# other functions found in thermal_plant_actor_ext.jl to be loaded if ThermalSystemModels is imported