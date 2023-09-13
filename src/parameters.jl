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

"""
    case_parameters_creation_from_ini(ini_case::AbstractParameters)

Useful function to create a case_paramter function for a case from an ini
"""
function case_parameters_creation_from_ini(ini_case::ParametersAllInits)
    ini_base = ParametersInits()
    for field in keys(ini_base)
        code_string = ""
        for sub_field in keys(getproperty(ini_base, field))
            value_case = getproperty(getproperty(ini_case, field), sub_field, missing)
            value_base = getproperty(getproperty(ini_base, field), sub_field, missing)
            if value_case !== value_base
                if typeof(value_case) <: Symbol
                    code_string = "ini.$field.$sub_field = :$value_case"
                elseif typeof(value_case) <: String
                    code_string = "ini.$field.$sub_field = \"$value_case\""
                elseif typeof(value_case) <: OrderedCollections.OrderedDict
                    code_string = "ini.$field.$sub_field = OrderedCollections.OrderedDict(\n    " * join([":$k => $v," for (k, v) in value_case], "\n    ") * "\n)"
                else
                    code_string = "ini.$field.$sub_field = $value_case"
                end
                println(code_string)
            end
        end
        if !isempty(code_string)
            println()
        end
    end
end