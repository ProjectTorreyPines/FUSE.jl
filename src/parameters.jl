using SimulationParameters

function SimulationParameters.Entry{T}(ids::Type, field::Symbol; default=missing, check=nothing) where {T}
    nfo = IMAS.info(ids, field)
    return Entry{T}(nfo.units, nfo.documentation; default, check)
end

function SimulationParameters.Switch{T}(options, ids::Type{<:IMAS.IDS}, field::Symbol; default=missing) where {T}
    nfo = IMAS.info(ids, field)
    return Switch{T}(options, nfo.units, nfo.documentation; default)
end

abstract type ParametersActor{T} <: AbstractParameters{T} end # abstract type for parameters of an actor (each actor has its concrete type)
abstract type ParametersAllActors{T} <: AbstractParameters{T} end # --> abstract type for all parameters of all actors (the concrete type of this is ParametersActors)

abstract type ParametersInit{T} <: AbstractParameters{T} end # abstract type for parameters of a init (each init has its concrete type)
abstract type ParametersAllInits{T} <: AbstractParameters{T} end # --> abstract type for all parameters of all inits (the concrete type of this is ParametersInits)

function SimulationParameters.time_range(parameters::ParametersInit)
    return SimulationParameters.time_range(SimulationParameters.top(parameters))
end

function SimulationParameters.time_range(parameters::ParametersAllInits)
    if ismissing(parameters.time, :pulse_shedule_time_basis)
        return Float64[]
    else
        return parameters.time.pulse_shedule_time_basis
    end
end

function SimulationParameters.global_time(parameters::ParametersInit)
    return SimulationParameters.global_time(SimulationParameters.top(parameters))
end

function SimulationParameters.global_time(parameters::ParametersAllInits)
    return parameters.time.simulation_start
end

function SimulationParameters.global_time(parameters::ParametersInit, time0::Float64)
    return SimulationParameters.global_time(SimulationParameters.top(parameters), time0)
end

function SimulationParameters.global_time(parameters::ParametersAllInits, time0::Float64)
    return parameters.time.simulation_start = time0
end

function Base.show(io::IO, ::MIME"text/plain", x::Tuple{Vararg{AbstractParameters}})
    println(io, join(map(SimulationParameters.spath,x),", "))
end
