using InteractiveUtils: subtypes
import AbstractTrees

abstract type Parameter end
abstract type Parameters end

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

function Switch(options::Vector{<:Pair}, units::String, description::String; default=missing)
    opts = Dict{Any,SwitchOption}()
    for (key, desc) in options
        opts[key] = SwitchOption(key, desc)
    end
    return Switch(missing, WeakRef(missing), opts, units, description, default, default, default)
end

function Switch(options::Vector{<:Union{Symbol,String}}, units::String, description::String; default=missing)
    opts = Dict{eltype(options),SwitchOption}()
    for key in options
        opts[key] = SwitchOption(key, "$key")
    end
    return Switch(missing, WeakRef(missing), opts, units, description, default, default, default)
end

function Switch(options, ids::Type{<:IMAS.IDS}, field::Symbol; default=missing)
    location = "$(IMAS.fs2u(ids)).$(field)"
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

function parameter_color(p::Parameter)
    value = p.value
    if value === missing
        color = :yellow
    elseif typeof(value) == typeof(p.default) && value == p.default
        color = :green
    elseif typeof(value) == typeof(p.base) && value == p.base
        color = :blue
    else
        color = :red
    end
end

function Base.show(io::IO, p::Parameter)
    color = parameter_color(p)
    printstyled(io, join(path(p), "."); bold=true, color=color)
    for item in fieldnames(typeof(p))
        if startswith(string(item), "_")
            continue
        end
        printstyled(io, "\n- $item: "; bold=true)
        printstyled(io, "$(getfield(p, item))")
    end
end

#= ====================== =#
#  Optimization parameter  #
#= ====================== =#
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

#= ============== =#
#  ParametersInit  #
#= ============== =#
mutable struct ParametersInit <: Parameters
    _name::Union{Missing,Symbol}
    _parent::WeakRef
    _parameters::Dict{Symbol,Union{Parameter,Parameters}}
end

function ParametersInit(::Nothing)
    return ParametersInit(missing, WeakRef(missing), Dict{Symbol,Union{Parameter,ParametersInit}}())
end

function ParametersInit(group::Symbol; kw...)
    if length(methods(ParametersInit, (Type{Val{group}},))) == 0
        throw(InexistentParameterException(ParametersInit, [group]))
    end
    par = ParametersInit(Val{group}; kw...)
    par._name = group
    return par
end

#= =============== =#
#  ParametersActor  #
#= =============== =#
mutable struct ParametersActor <: Parameters
    _name::Union{Missing,Symbol}
    _parent::WeakRef
    _parameters::Dict{Symbol,Union{Parameter,Parameters}}
end

function ParametersActor(::Nothing)
    return ParametersActor(missing, WeakRef(missing), Dict{Symbol,Union{Parameter,ParametersActor}}())
end

function ParametersActor(group::Symbol; kw...)
    if length(methods(ParametersActor, (Type{Val{group}},))) == 0
        throw(InexistentParameterException(ParametersActor, [group]))
    end
    pars = ParametersActor(Val{group}; kw...)
    pars._name = group
    return pars
end

"""
    ParametersActor()

Generates actor parameters 
"""
function ParametersActor()
    act = ParametersActor(missing, WeakRef(missing), Dict{Symbol,Union{Parameter,ParametersActor}}())
    for par in concretetypes(AbstractActor)
        par = Symbol(replace(string(par), "FUSE." => ""))
        try
            setproperty!(act, par, ParametersActor(par))
        catch e
            if typeof(e) <: InexistentParameterException
                @warn e
            else
                rethrow()
            end
        end
    end

    act._name = :act
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

"""
    getproperty(p::Parameters, key::Symbol, default)

Return value of `key` parameter or `default` if parameter is missing
NOTE: This is useful because accessing a `missing` parameter would raise an error
"""
function Base.getproperty(p::Parameters, key::Symbol, default)
    value = p[key].value
    if value === missing
        return default
    else
        return value
    end
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

function Base.iterate(par::FUSE.Parameters)
    Base.iterate(par, collect(keys(par)))
end

function Base.iterate(par::FUSE.Parameters, state)
    if isempty(state)
        return nothing
    end
    key = popfirst!(state)
    data = par[key].value
    return key => data, state
end

function Base.show(io::IO, ::MIME"text/plain", pars::Parameters, depth::Int=0)
    return AbstractTrees.print_tree(io, pars)
end

function AbstractTrees.children(pars::Parameters)
    return [pars[k] for k in sort(collect(keys(pars)))]
end

function AbstractTrees.printnode(io::IO, par::Parameter)
    color = parameter_color(par)
    printstyled(io, par._name)
    printstyled(io, " ➡ ")
    printstyled(io, "$(repr(par.value))"; color=color)
    if length(par.units) > 0 && par.value !== missing
        printstyled(io, " [$(par.units)]"; color=color)
    end
end

function AbstractTrees.printnode(io::IO, pars::Parameters)
    printstyled(io, pars._name; bold=true)
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

function doc(parameters::Parameters)
    if typeof(parameters) <: ParametersActor
        ppath = "act.$(parameters._name)"
    else
        ppath = "ini.$(parameters._name)"
    end
    txt = []
    for par in sort(collect(keys(parameters)))
        if typeof(parameters[par]) <: Parameters
            push!(txt, "**`$(ppath).$par`**: $(typeof(parameters[par]))")
        else
            if isempty(parameters[par].units)
                units = ""
            else
                units = " [$(parameters[par].units)]"
            end
            push!(txt, "**`$(ppath).$par`**:$units $(parameters[par].description)")
        end
    end
    if isempty(txt)
        return ""
    else
        return "* " * join(txt, "\n* ")
    end
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
Base.showerror(io::IO, e::BadParameterException) =
    print(io, "ERROR: Parameter $(join(e.path,".")) = $(repr(e.value)) is not one of the valid options: $(join(map(repr,e.options),", "))")

#= ============ =#
#  case studies  #
#= ============ =#
# NOTE only called once at precompile time, kernel needs to be restarted to include new file in cases
for filename in readdir(joinpath(dirname(@__FILE__), "..", "cases"))
    if endswith(filename, ".jl")
        include("../cases/" * filename)
    end
end

function case_parameters(case::Symbol; kw...)
    if length(methods(case_parameters, (Type{Val{case}},))) == 0
        throw(InexistentParameterException(Parameters, [case]))
    end
    return case_parameters(Val{case}; kw...)
end