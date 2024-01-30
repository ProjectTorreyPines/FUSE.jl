using SimulationParameters

function SimulationParameters.Entry{T}(ids::Type, field::Symbol; default=missing, check=nothing) where {T}
    txt = IMAS.info(ids, field)
    return Entry{T}(get(txt, "units", "-"), get(txt, "documentation", ""); default, check)
end

function SimulationParameters.Switch{T}(options, ids::Type{<:IMAS.IDS}, field::Symbol; default=missing) where {T}
    txt = IMAS.info(ids, field)
    return Switch{T}(options, get(txt, "units", "-"), get(txt, "documentation", ""); default)
end

abstract type ParametersActor <: AbstractParameters end # container for all parameters of an actor
abstract type ParametersAllActors <: AbstractParameters end # --> abstract type of ParametersActors, container for all parameters of all actors

abstract type ParametersInit <: AbstractParameters end # container for all parameters of a init
abstract type ParametersAllInits <: AbstractParameters end # --> abstract type of ParametersInits, container for all parameters of all inits
