abstract type Parameter end

#= ================= =#
#   Scalar Parameter  #
#= ================= =#
mutable struct Entry{T} <: Parameter
    units::String
    description::String
    value::T
    base::T
    default::T
end

"""
    Entry(T, units::String, description::String; default = missing)

Defines a parameter
"""
function Entry(T, units::String, description::String; default = missing)
    return Entry{Union{Missing,T}}(units, description, default, default, default)
end

function Entry(T, ids::Type, field::Symbol; default = missing)
    info = IMAS.imas_info(ids, field)
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
    units::String
    description::String
    value::Union{Missing,Symbol}
    base::Union{Missing,Symbol}
    default::Union{Missing,Symbol}
end

"""
    Switch(options, units::String, description::String; default = missing)

Defines a switch
"""
function Switch(options::Dict{Symbol,Union{SwitchOption,Entry}}, units::String, description::String; default = missing)
    if !in(default, keys(options))
        error("$(repr(default)) is not a valid option: $(collect(keys(options)))")
    end
    return Switch(options, units, description, default, default, default)
end

function Switch(options::Vector{T}, units::String, description::String; default = missing) where {T<:Pair{Symbol,Z}} where {Z<:Any}
    opts = Dict{Symbol,SwitchOption}()
    for (key, desc) in options
        opts[key] = SwitchOption(key, units, desc)
    end
    return Switch(opts, units, description, default, default, default)
end

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

struct BadParameterException <: Exception
    key::Symbol
    valid::Vector{Symbol}
end
Base.showerror(io::IO, e::BadParameterException) = print(io, "ERROR: `$(repr(e.key))` is not one of the valid options: $(e.valid)")

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
        if (value !== missing) && !(value in keys(parameter.options))
            throw(BadParameterException(value, collect(keys(parameter.options))))
        end
        return parameter.value = value
    else
        error("Unrecognized type $(typeof(parameter))")
    end
end

function Base.show(io::IO, p::Parameters, depth::Int)
    _parameters = getfield(p, :_parameters)
    for item in sort(collect(keys(_parameters)))
        parameter = _parameters[item]
        if typeof(parameter) <: Parameters
            printstyled(io, "$(' '^(2*depth))")
            printstyled(io, "$(item)\n"; bold = true)
            show(io, parameter, depth + 1)
        else
            value = parameter.value
            units = parameter.units
            if value === missing
                continue
            elseif value === parameter.default
                color = :blue
            elseif value === parameter.base
                color = :black
            else
                color = :red
            end
            printstyled(io, "$(' '^(2*depth))")
            printstyled(io, "$(item)"; color = color)
            printstyled(io, " ➡ "; color = :red)
            printstyled(io, "$(value)"; color = color)
            if length(units) > 0
                printstyled(io, " [$(units)]"; color = color)
            end
            printstyled(io, "\n")
        end
    end
end

function set_new_base(p::Parameters)
    _parameters = getfield(p, :_parameters)
    for item in sort(collect(keys(_parameters)))
        parameter = _parameters[item]
        if typeof(parameter) <: Parameters
            set_new_base(parameter)
        else
            parameter.base = parameter.value
        end
    end
    return p
end

function Base.show(io::IO, ::MIME"text/plain", p::Parameters)
    return show(io, p, 0)
end

function Base.ismissing(p::Parameters, field::Symbol)::Bool
    return getfield(p, :_parameters)[field] === missing
end

