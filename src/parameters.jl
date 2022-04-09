abstract type Parameter end

#= ===== =#
#  Entry  #
#= ===== =#
mutable struct Entry{T} <: Parameter
    units::String
    description::String
    value::T
    base::T
    default::T
end

"""
    Entry(T, units::String, description::String; default = missing)

Defines a entry parameter
"""
function Entry(T, units::String, description::String; default = missing)
    return Entry{Union{Missing,T}}(units, description, default, default, default)
end

function Entry(T, ids::Type, field::Symbol; default = missing)
    txt = IMAS.info(ids, field)
    return Entry(T, get(txt, "units", ""), get(txt, "documentation", ""); default)
end

#= ====== =#
#  Switch  #
#= ====== =#

struct SwitchOption
    value::Any
    description::String
end

mutable struct Switch <: Parameter
    options::Dict{T,SwitchOption} where {T<:Any}
    units::String
    description::String
    value::Any
    base::Any
    default::Any
end

"""
    Switch(options, units::String, description::String; default = missing)

Defines a switch parameter
"""
function Switch(options::Dict{Any,SwitchOption}, units::String, description::String; default = missing)
    if !in(default, keys(options))
        error("$(repr(default)) is not a valid option: $(collect(keys(options)))")
    end
    return Switch(options, units, description, default, default, default)
end

function Switch(options::Vector{T}, units::String, description::String; default = missing) where {T<:Pair}
    opts = Dict{Any,SwitchOption}()
    for (key, desc) in options
        opts[key] = SwitchOption(key, desc)
    end
    return Switch(opts, units, description, default, default, default)
end

function Switch(options::Vector{T}, units::String, description::String; default = missing) where {T<:Union{Symbol,String}}
    opts = Dict{eltype(options),SwitchOption}()
    for key in options
        opts[key] = SwitchOption(key, "$key")
    end
    return Switch(opts, units, description, default, default, default)
end

function Switch(options, ids::Type{T}, field::Symbol; default = missing) where {T<:IMAS.IDS}
    location = "$(IMAS._f2u(ids)).$(field)"
    txt = IMAS.info(location)
    return Switch(options, get(txt, "units", ""), get(txt, "documentation", ""); default)
end

#= ========== =#
#  Parameters  #
#= ========== =#

abstract type Parameters end

mutable struct InitParameters <: Parameters
    _path::Vector{Symbol}
    _parameters::Dict{Symbol,Union{Parameter,Parameters}}
end

function InitParameters(::Nothing)
    return InitParameters(Symbol[], Dict{Symbol,Union{Parameter,InitParameters}}())
end

function InitParameters(group::Symbol; kw...)
    if length(methods(InitParameters, (Type{Val{group}},))) == 0
        throw(InexistentParameterException(InitParameters, [group]))
    end
    return InitParameters(Val{group}; kw...)
end

mutable struct ModelParameters <: Parameters
    _path::Vector{Symbol}
    _parameters::Dict{Symbol,Union{Parameter,Parameters}}
end

function ModelParameters(::Nothing)
    return ModelParameters(Symbol[], Dict{Symbol,Union{Parameter,ModelParameters}}())
end

function ModelParameters(group::Symbol; kw...)
    if length(methods(ModelParameters, (Type{Val{group}},))) == 0
        throw(InexistentParameterException(ModelParameters, [group]))
    end
    return ModelParameters(Val{group}; kw...)
end

function Base.fieldnames(p::Parameters)
    return collect(keys(getfield(p, :_parameters)))
end

function Base.getproperty(p::Parameters, key::Symbol)
    _parameter = getfield(p, :_parameters)
    if !(key in keys(_parameter))
        throw(InexistentParameterException(typeof(p), vcat(getfield(p, :_path), key)))
    end
    parameter = _parameter[key]

    if typeof(parameter) <: Parameters
        value = parameter
    elseif typeof(parameter) <: Entry
        value = parameter.value
    elseif typeof(parameter) <: Switch
        if parameter.value === missing
            throw(NotsetParameterException(vcat(getfield(p, :_path), key), collect(keys(parameter.options))))
        end
        value = parameter.options[parameter.value].value
    else
        error("Unrecognized type $(typeof(parameter))")
    end

    if value === missing
        throw(NotsetParameterException(vcat(getfield(p, :_path), key)))
    end

    return value
end

function Base.setproperty!(p::Parameters, key::Symbol, value)
    if typeof(value) <: Union{Parameter,Parameters}
        if typeof(value) <: Parameters
            setfield!(value, :_path, vcat(getfield(p, :_path), key))
        end
        getfield(p, :_parameters)[key] = value
        return value
    end

    _parameter = getfield(p, :_parameters)
    if !(key in keys(_parameter))
        throw(InexistentParameterException(typeof(p), vcat(getfield(p, :_path), key)))
    end
    parameter = _parameter[key]

    if typeof(parameter) <: Switch
        try
            return parameter.value = value
        catch e
            if typeof(e) <: BadParameterException # retrhow the exception but add more to the path information
                throw(BadParameterException(vcat(getfield(p, :_path), key), value, collect(keys(parameter.options))))
            end
        end
    else
        return parameter.value = value
    end

    return value
end

function Base.setproperty!(p::Switch, key::Symbol, value)
    if typeof(value) <: Pair
        p.options[value.first].value = value.second
        value = value.first
    end
    if (value !== missing) && !(value in keys(p.options))
        throw(BadParameterException([key], value, collect(keys(p.options))))
    end
    return setfield!(p, :value, value)
end

function Base.show(io::IO, p::Parameters, depth::Int)
    _parameters = getfield(p, :_parameters)
    for item in sort(collect(keys(_parameters)))
        parameter = _parameters[item]
        if typeof(parameter) <: Parameters
            printstyled(io, "$(" "^(2*depth))")
            printstyled(io, "$(item)\n"; bold = true)
            show(io, parameter, depth + 1)
        else
            value = parameter.value
            units = parameter.units
            if value === missing
                color = :yellow
            elseif typeof(value) == typeof(parameter.default) && value == parameter.default
                color = :green
            elseif typeof(value) == typeof(parameter.base) && value == parameter.base
                color = :blue
            else
                color = :red
            end
            printstyled(io, "$(" "^(2*depth))")
            printstyled(io, "$(item)"; color = color)
            printstyled(io, " ➡ "; color = :red)
            printstyled(io, "$(repr(value))"; color = color)
            if length(units) > 0 && value !== missing
                printstyled(io, " [$(units)]"; color = color)
            end
            printstyled(io, "\n")
        end
    end
end

function set_new_base!(p::Parameters)
    _parameters = getfield(p, :_parameters)
    for item in sort(collect(keys(_parameters)))
        parameter = _parameters[item]
        if typeof(parameter) <: Parameters
            set_new_base!(parameter)
        else
            setfield!(parameter, :base, parameter.value)
        end
    end
    return p
end

function Base.show(io::IO, ::MIME"text/plain", p::Parameters)
    return show(io, p, 0)
end

function Base.ismissing(p::Parameters, field::Symbol)::Bool
    return getfield(p, :_parameters)[field].value === missing
end

#= ================= =#
#  Parameters errors  #
#= ================= =#
struct InexistentParameterException <: Exception
    parameter_type::DataType
    path::Vector{Symbol}
end
Base.showerror(io::IO, e::InexistentParameterException) = print(io, "ERROR: $(e.parameter_type) $(join(e.path,".")) does not exist")

struct NotsetParameterException <: Exception
    path::Vector{Symbol}
    options::Vector{Any}
end
NotsetParameterException(path::Vector{Symbol}) = NotsetParameterException(path, [])
function Base.showerror(io::IO, e::NotsetParameterException)
    if length(e.options) > 0
        print(io, "ERROR: Parameter $(join(e.path,".")) is not set. Valid options are: $(join(map(repr,e.options),", "))")
    else
        print(io, "ERROR: Parameter $(join(e.path,".")) is not set")
    end
end

struct BadParameterException <: Exception
    path::Vector{Symbol}
    value::Any
    options::Vector{Any}
end
Base.showerror(io::IO, e::BadParameterException) = print(io, "ERROR: Parameter $(join(e.path,".")) value `$(repr(e.value))` is not one of the valid options: $(join(map(repr,e.options),", "))")
