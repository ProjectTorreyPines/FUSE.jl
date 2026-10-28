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

"""
    @actor_parameters_struct ActorName{T} body

Generate a mutable struct for actor parameters with standard boilerplate fields automatically included.

The macro creates a struct named `FUSEparameters__ActorName` (or `_FUSEparameters__ActorName` for underscore-prefixed actors)
that inherits from `ParametersActor{T}` and includes the standard fields:
- `_parent::WeakRef = WeakRef(nothing)`
- `_name::Symbol = :not_set` 
- `_time::Float64 = NaN`
"""
macro actor_parameters_struct(name_expr, body)
    # Parse the actor name and type parameter from expressions like ActorName{T}
    if name_expr isa Expr && name_expr.head == :curly
        actor_name = name_expr.args[1]
        type_param = name_expr.args[2]
    else
        error("Expected format: ActorName{T}")
    end
    
    # Convert actor name to string and determine struct name
    actor_str = string(actor_name)
    if startswith(actor_str, "_")
        # Remove _ in underscore-prefixed actors
        struct_name = Symbol("_FUSEparameters__" * actor_str[2:end])
    else
        struct_name = Symbol("FUSEparameters__" * actor_str)
    end
    
    # Extract field definitions from the body
    user_fields = if body isa Expr && body.head == :block
        filter(x -> !(x isa LineNumberNode), body.args)
    else
        [body]
    end
    
    # Standard boilerplate fields that every actor needs
    standard_fields = [
        :(_parent::WeakRef = WeakRef(nothing)),
        :(_name::Symbol = :not_set),
        :(_time::Float64 = NaN)
    ]
    
    # Combine standard fields with user-defined fields
    all_fields = vcat(standard_fields, user_fields)

    # Resolve the supertype here (in this module) and splice the type object
    # into the expansion. Otherwise the bare symbol `ParametersActor` inside
    # the escaped quote would be looked up in the caller module's scope,
    # forcing every caller to `using FUSE: ParametersActor`.
    par_actor_type = ParametersActor

    # Generate the @kwdef struct
    return esc(quote
        Base.@kwdef mutable struct $struct_name{$type_param<:Real} <: $par_actor_type{$type_param}
            $(all_fields...)
        end
    end)
end

"""
    unset_parameters!(par::AbstractParameters)

Set every parameter of `par` to `missing` (parameter vectors are not visited)
"""
function unset_parameters!(par::AbstractParameters)
    for field in keys(par)
        p = SimulationParameters.getparameter(par, field)
        if p isa AbstractParameters
            unset_parameters!(p)
        elseif p isa SimulationParameters.AbstractParameter
            setfield!(p, :value, missing)
        end
    end
    return par
end

"""
    apply_act_overrides!(act::ParametersAllActors, overrides::ParametersAllActors)

Set on `act` every actor parameter that is set (not `missing`) in `overrides`, e.g. one started with [`unset_parameters!`](@ref)
"""
function apply_act_overrides!(act::ParametersAllActors, overrides::ParametersAllActors)
    _apply_overrides!(act, overrides)
    return act
end

function _apply_overrides!(par::AbstractParameters, overrides::AbstractParameters)
    for field in keys(overrides)
        p = SimulationParameters.getparameter(overrides, field)
        if p isa AbstractParameters
            _apply_overrides!(getproperty(par, field), p)
        elseif p isa SimulationParameters.AbstractParameter && !ismissing(p.value)
            setproperty!(par, field, p.value)
        end
    end
end
