abstract type Parameter end

#= ================= =#
#   Scalar Parameter  #
#= ================= =#
mutable struct Entry{T} <: Parameter
    default::T
    units::String
    description::String
    value::T
end

"""
    Entry(default, units::String, description::String)

Defines a parameter
"""
function Entry(T, units::String, description::String; default = missing)
    return Entry{Union{Missing,T}}(default, units, description, default)
end

"""
    Entry(default, ids::Type{T}, field::Symbol) where {T<:IMAS.IDS}

Defines a parameter taking units and description from IDS field
"""
function Entry(T, ids, field::Symbol; default = missing)
    location = "$(IMAS._f2u(ids)).$(field)"
    info = IMAS.imas_info(location)
    return Entry(T, get(info, "units", ""), get(info, "documentation", ""); default)
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
#   Switch   #
#= ================= =#

mutable struct Switch <: Parameter
    options::Dict{Symbol,Union{SwitchOption,Entry}}
    default::Symbol
    units::String
    description::String
    value::Symbol
end
"""
    Switch(options, default, units::String, description::String)

Defines a switch
"""
function Switch(options, units::String, description::String; default = missing)
    if !in(default, keys(options))
        error("$(repr(default)) is not a valid option: $(collect(keys(options)))")
    end
    return Switch(options, default, units, description, default)
end

"""
    Switch(options, default, ids::IMAS.IDS, field::Symbol)

Defines a switch taking units and description from IDS field
"""
function Switch(options, ids::Type{T}, field::Symbol; default = missing) where {T<:IMAS.IDS}
    location = "$(IMAS._f2u(ids)).$(field)"
    info = IMAS.imas_info(location)
    return Switch(options, get(info, "units", ""), get(info, "documentation", ""); default)
end

#= ================= =#
#   Fuse Parameters   #
#= ================= =#

struct InexistentParameterException <: Exception
    key::Symbol
end
Base.showerror(io::IO, e::InexistentParameterException) = print(io, "ERROR: parameter $(e.key) does not exist")

struct NotsetParameterException <: Exception
    key::Symbol
end
Base.showerror(io::IO, e::NotsetParameterException) = print(io, "ERROR: parameter $(e.key) is not set")

struct Parameters
    _parameters::Dict{Symbol,Union{Parameter,Parameters}}
end

function Base.fieldnames(p::Parameters)
    return collect(keys(getfield(p, :_parameters)))
end


function Base.getproperty(p::Parameters, key::Symbol)
    _parameter = getfield(p, :_parameters)
    if !(key in keys(_parameter))
        throw(InexistentParameterException(key))
    end
    parameter = _parameter[key]

    if typeof(parameter) <: Parameters
        value = parameter
    elseif typeof(parameter) <: Entry
        value = parameter.value
    elseif typeof(parameter) <: Switch
        value = parameter.options[parameter.value].value
    else
        error("Unrecognized type $(typeof(parameter))")
    end

    if ismissing(value)
        throw(NotsetParameterException(key))
    end

    return value
end


function Base.setproperty!(p::Parameters, key::Symbol, value)
    if typeof(value) <: Union{Parameter,Parameters}
        getfield(p, :_parameters)[key] = value
        return
    end

    _parameter = getfield(p, :_parameters)
    if !(key in keys(_parameter))
        throw(InexistentParameterException(key))
    end
    parameter = _parameter[key]

    if typeof(parameter) <: Entry
        return parameter.value = value
    elseif typeof(parameter) <: Switch
        if typeof(value) <: Pair
            parameter.options[value.first].value = value.second
            value = value.first
        end
        if !in(value, keys(parameter.options))
            error("$(repr(value)) is not a valid option: $(collect(keys(parameter.options)))")
        end
        return parameter.value = value
    else
        error("Unrecognized type $(typeof(parameter))")
    end
end


#= ================= =#
#   Parameters list   #
#= ================= =#

include("parameters_list.jl")
