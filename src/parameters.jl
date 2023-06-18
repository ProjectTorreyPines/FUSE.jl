using SimulationParameters

function SimulationParameters.Entry{T}(ids::Type, field::Symbol; default=missing) where {T}
    txt = IMAS.info(ids, field)
    return Entry{T}(get(txt, "units", "-"), get(txt, "documentation", ""); default)
end

function SimulationParameters.Switch{T}(options, ids::Type{<:IMAS.IDS}, field::Symbol; default=missing) where {T}
    txt = IMAS.info(ids, field)
    return Switch{T}(options, get(txt, "units", "-"), get(txt, "documentation", ""); default)
end

abstract type ParametersActor <: AbstractParameters end # container for all parameters of an actor
abstract type ParametersAllActors <: AbstractParameters end # --> abstract type of ParametersActors, container for all parameters of all actors

abstract type ParametersInit <: AbstractParameters end # container for all parameters of a init
abstract type ParametersAllInits <: AbstractParameters end # --> abstract type of ParametersInits, container for all parameters of all inits

#= ============ =#
#  case studies  #
#= ============ =#
# NOTE only called once at precompile time, kernel needs to be restarted to include new file in cases
for filename in readdir(joinpath(@__DIR__, "..", "cases"))
    if endswith(filename, ".jl")
        include("../cases/" * filename)
    end
end

function case_parameters(case::Symbol; kw...)
    if length(methods(case_parameters, (Type{Val{case}},))) == 0
        throw(InexistentParameterException([case]))
    end
    return case_parameters(Val{case}; kw...)
end
