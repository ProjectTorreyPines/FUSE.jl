abstract type AbstractParameter end

#= ================= =#
#   Scalar Parameter  #
#= ================= =#
mutable struct ScalarParameter{T} <: AbstractParameter
    default::T
    units::String
    description::String
    value::T
end

"""
    ScalarParameter(default, units, description)

Defines a scalar parameter
"""
function ScalarParameter(default, units, description)
    return ScalarParameter(default, units, description, default)
end

#= =============== =#
#   Switch Option   #
#= =============== =#

struct SwitchOption
    value::Any
    units::String
    description::String
end

#= ================= =#
#   SwitchParameter   #
#= ================= =#

mutable struct SwitchParameter <: AbstractParameter
    options::Dict{Symbol,Union{SwitchOption,ScalarParameter}}
    default::Symbol
    description::String
    value::Symbol
end


function SwitchParameter(options, default, description)
    if ! in(default, keys(options))
        error("$(repr(default)) is not a valid option: $(collect(keys(options)))")
    end
    return SwitchParameter(options, default, description, default)
end

#= ================= =#
#   Fuse Parameters   #
#= ================= =#

mutable struct FuseParameters
    parameters::Dict{Symbol,AbstractParameter}
end


function Base.getindex(p::FuseParameters, key)
    parameter = p.parameters[key]
    if typeof(parameter) <: ScalarParameter
        return parameter.value
    elseif typeof(parameter) <: SwitchParameter
        parameter.options[parameter.value].value
    else
        throw(KeyError(key))
    end
end


function Base.setindex!(p::FuseParameters, value, key)
    parameter = p.parameters[key]
    if typeof(parameter) <: ScalarParameter
        return parameter.value = value
    elseif typeof(parameter) <: SwitchParameter
        if typeof(value) <: Pair
            parameter.options[value.first].value = value.second
            value = value.first
        end
        if ! in(value, keys(parameter.options))
            error("$(repr(value)) is not a valid option: $(collect(keys(parameter.options)))")
        end
        return parameter.value = value
    else
        throw(KeyError(key))
    end
end


#= ================= =#
#   Parameters list   #
#= ================= =#

include("parameters_list.jl")
