using InteractiveUtils: subtypes

abstract type Parameter end
abstract type Parameters end

#= ======================= =#
#  Optimization parameters  #
#= ======================= =#
struct OptParameter
    nominal::Real
    lower::Real
    upper::Real
end

function ↔(x::Real, r::AbstractVector)
    @assert typeof(x) == typeof(r[1]) == typeof(r[end]) "type of optimization range does not match the nominal value"
    return OptParameter(x, r[1], r[end])
end

function opt_parameters(p::Parameters, opt_vector=Parameter[])
    _parameters = getfield(p, :_parameters)
    for k in keys(_parameters)
        parameter = _parameters[k]
        if typeof(parameter) <: Parameters
            opt_parameters(parameter, opt_vector)
        elseif typeof(parameter) <: Entry
            if parameter.lower !== missing
                push!(opt_vector, parameter)
            end
        end
    end
    return opt_vector
end

#= ===== =#
#  Entry  #
#= ===== =#
mutable struct Entry{T} <: Parameter
    _name::Union{Missing,Symbol}
    _parent::WeakRef
    units::String
    description::String
    value::T
    base::T
    default::T
    lower::Union{Missing,Float64}
    upper::Union{Missing,Float64}
end

"""
    Entry(T, units::String, description::String; default = missing)

Defines a entry parameter
"""
function Entry(T, units::String, description::String; default=missing)
    return Entry{Union{Missing,T}}(missing, WeakRef(missing), units, description, default, default, default, missing, missing)
end

function Entry(T, ids::Type, field::Symbol; default=missing)
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
    _name::Union{Missing,Symbol}
    _parent::WeakRef
    options::Dict{Any,SwitchOption}
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
function Switch(options::Dict{Any,SwitchOption}, units::String, description::String; default=missing)
    if !in(default, keys(options))
        error("$(repr(default)) is not a valid option: $(collect(keys(options)))")
    end
    return Switch(missing, WeakRef(missing), options, units, description, default, default, default)
end

function Switch(options::Vector{T}, units::String, description::String; default=missing) where {T<:Pair}
    opts = Dict{Any,SwitchOption}()
    for (key, desc) in options
        opts[key] = SwitchOption(key, desc)
    end
    return Switch(missing, WeakRef(missing), opts, units, description, default, default, default)
end

function Switch(options::Vector{T}, units::String, description::String; default=missing) where {T<:Union{Symbol,String}}
    opts = Dict{eltype(options),SwitchOption}()
    for key in options
        opts[key] = SwitchOption(key, "$key")
    end
    return Switch(missing, WeakRef(missing), opts, units, description, default, default, default)
end

function Switch(options, ids::Type{T}, field::Symbol; default=missing) where {T<:IMAS.IDS}
    location = "$(IMAS._f2u(ids)).$(field)"
    txt = IMAS.info(location)
    return Switch(options, get(txt, "units", ""), get(txt, "documentation", ""); default)
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

#= ============== =#
#  InitParameters  #
#= ============== =#
mutable struct InitParameters <: Parameters
    _name::Union{Missing,Symbol}
    _parent::WeakRef
    _parameters::Dict{Symbol,Union{Parameter,Parameters}}
end

function InitParameters(::Nothing)
    return InitParameters(missing, WeakRef(missing), Dict{Symbol,Union{Parameter,InitParameters}}())
end

function InitParameters(group::Symbol; kw...)
    if length(methods(InitParameters, (Type{Val{group}},))) == 0
        throw(InexistentParameterException(InitParameters, [group]))
    end
    return InitParameters(Val{group}; kw...)
end

#= =============== =#
#  ActorParameters  #
#= =============== =#
mutable struct ActorParameters <: Parameters
    _name::Union{Missing,Symbol}
    _parent::WeakRef
    _parameters::Dict{Symbol,Union{Parameter,Parameters}}
end

function ActorParameters(::Nothing)
    return ActorParameters(missing, WeakRef(missing), Dict{Symbol,Union{Parameter,ActorParameters}}())
end

function ActorParameters(group::Symbol; kw...)
    if length(methods(ActorParameters, (Type{Val{group}},))) == 0
        throw(InexistentParameterException(ActorParameters, [group]))
    end
    return ActorParameters(Val{group}; kw...)
end

"""
    ActorParameters()

Generates actor parameters 
"""
function ActorParameters()
    act = ActorParameters(missing, WeakRef(missing), Dict{Symbol,Union{Parameter,ActorParameters}}())
    for par in subtypes(ActorAbstract)
        par = Symbol(replace(string(par), "FUSE." => ""))
        try
            setproperty!(act, par, ActorParameters(par))
        catch e
            if typeof(e) <: InexistentParameterException
                @warn e
            else
                rethrow()
            end
        end
    end
    return act
end

#= ========== =#
#  Parameters  #
#= ========== =#
function path(p::Union{Parameter,Parameters})
    name = getfield(p, :_name)
    if name === missing
        return Symbol[]
    end
    pp = Symbol[name]
    while typeof(p._parent.value) <: Parameters
        if p._parent.value._name === missing
            break
        end
        pushfirst!(pp, p._parent.value._name)
        p = p._parent.value
    end
    return pp
end

function Base.keys(p::Parameters)
    return keys(getfield(p, :_parameters))
end

function Base.getindex(p::Parameters, field::Symbol)
    return getfield(p, :_parameters)[field]
end

function Base.setindex!(p::Parameters, value::Any, field::Symbol)
    return getfield(p, :_parameters)[field] = value
end

function Base.getproperty(p::Parameters, key::Symbol)
    if key ∈ fieldnames(typeof(p))
        return getfield(p, key)
    elseif key ∉ keys(p)
        throw(InexistentParameterException(typeof(p), vcat(path(p), key)))
    end
    parameter = p[key]

    if typeof(parameter) <: Parameters
        value = parameter
    elseif typeof(parameter) <: Entry
        value = parameter.value
    elseif typeof(parameter) <: Switch
        if parameter.value === missing
            throw(NotsetParameterException(vcat(path(p), key), collect(keys(parameter.options))))
        end
        value = parameter.options[parameter.value].value
    else
        error("Unrecognized type $(typeof(parameter))")
    end

    if value === missing
        throw(NotsetParameterException(vcat(path(p), key)))
    end

    return value
end

function Base.deepcopy(p::Union{Parameter,Parameters})
    p1 = Base.deepcopy_internal(p, Base.IdDict())
    p1._parent = WeakRef(missing)
    return p1
end

function Base.setproperty!(p::Parameters, key::Symbol, value)
    if key ∈ fieldnames(typeof(p))
        return setfield!(p, key, value)
    elseif typeof(value) <: Union{Parameter,Parameters}
        if typeof(value._parent.value) <: Union{Parameter,Parameters}
            value = deepcopy(value)
        end
        setfield!(value, :_parent, WeakRef(p))
        setfield!(value, :_name, key)
        p[key] = value
        return value
    end

    if !(key in keys(p))
        throw(InexistentParameterException(typeof(p), vcat(path(p), key)))
    end
    parameter = p[key]

    if typeof(parameter) <: Switch
        try
            return parameter.value = value
        catch e
            if typeof(e) <: BadParameterException # retrhow the exception but add more to the path information
                throw(BadParameterException(vcat(path(p), key), value, collect(keys(parameter.options))))
            end
        end
    else
        if typeof(value) <: OptParameter
            parameter.value = value.nominal
            if typeof(value.nominal) <: Integer
                parameter.lower = value.lower - 0.5
                parameter.upper = value.upper + 0.5
            else
                parameter.lower = value.lower
                parameter.upper = value.upper
            end
        else
            return parameter.value = value
        end
    end

    return value
end

function Base.show(io::IO, p::Parameters, depth::Int)
    for item in sort(collect(keys(p)))
        parameter = p[item]
        if typeof(parameter) <: Parameters
            printstyled(io, "$(" "^(2*depth))")
            printstyled(io, "$(item)\n"; bold=true)
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
            printstyled(io, "$(item)"; color=color)
            printstyled(io, " ➡ "; color=:red)
            printstyled(io, "$(repr(value))"; color=color)
            if length(units) > 0 && value !== missing
                printstyled(io, " [$(units)]"; color=color)
            end
            printstyled(io, "\n")
        end
    end
end

function set_new_base!(p::Parameters)
    for item in keys(p)
        parameter = p[item]
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
    return p[field].value === missing
end

"""
    (par::Parameters)(kw...)

This functor is used to override the parameters at function call
"""
function (par::Parameters)(kw...)
    if !isempty(kw)
        par = deepcopy(par)
        for (key, value) in kw
            setproperty!(par, key, value)
        end
    end
    return par
end

"""
    evalmissing(p::Parameters, field::Symbol) 

Return parameter value or `missing` if parameter is missing
NOTE: This is useful because accessing a `missing` parameter would raise an error
"""
function evalmissing(p::Parameters, field::Symbol)
    return p[field].value
end

"""
    par2dict(par::Parameters)

Convert FUSE parameters to dictionary
"""
function par2dict(par::Parameters)
    ret = Dict()
    return par2dict(par, ret)
end

function par2dict(par::Parameters, ret::AbstractDict)
    data = getfield(par, :_parameters)
    return par2dict(data, ret)
end

function par2dict(data::AbstractDict, ret::AbstractDict)
    for item in keys(data)
        if typeof(data[item]) <: Parameters
            ret[item] = Dict()
            par2dict(data[item], ret[item])
        elseif typeof(data[item]) <: Parameter
            ret[item] = Dict()
            ret[item][:value] = data[item].value
            ret[item][:units] = data[item].units
            ret[item][:description] = data[item].description
        end
    end
    return ret
end

"""
    par2json(@nospecialize(par::Parameters), filename::String; kw...)

Save the FUSE parameters to a JSON file with give `filename`
`kw` arguments are passed to the JSON.print function
"""
function par2json(@nospecialize(par::Parameters), filename::String; kw...)
    open(filename, "w") do io
        JSON.print(io, par2dict(par); kw...)
    end
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
Base.showerror(io::IO, e::BadParameterException) = print(io, "ERROR: Parameter $(join(e.path,".")) = $(repr(e.value)) is not one of the valid options: $(join(map(repr,e.options),", "))")

#= ============ =#
#  case studies  #
#= ============ =#
# NOTE only called once at precompile time, kernel needs to be restarted to include new file in cases
for filename in readdir(joinpath(dirname(@__FILE__), "..", "cases"))
    include("../cases/" * filename)
end

function case_parameters(case::Symbol; kw...)
    if length(methods(case_parameters, (Type{Val{case}},))) == 0
        throw(InexistentParameterException(case_parameters, [case]))
    end
    return case_parameters(Val{case}; kw...)
end