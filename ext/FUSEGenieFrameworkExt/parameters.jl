# Parameter introspection and manipulation for FUSE actors

"""
    ParamInfo

Information about a single parameter for UI rendering.

# Fields
- `name::Symbol` - Parameter name
- `param_type::Symbol` - :entry or :switch
- `value_type::String` - Type as string ("Int", "Float64", "Vector", "Symbol", "Bool", etc.)
- `current_value::Any` - Current value
- `default_value::Any` - Default value
- `options::Vector{Any}` - Available options (for switches only)
- `units::String` - Units string
- `description::String` - Parameter description
"""
struct ParamInfo
    name::Symbol
    param_type::Symbol  # :entry or :switch
    value_type::String  # "Int", "Float64", "Vector", "Symbol", "Bool"
    current_value::Any
    default_value::Any
    options::Vector{Any}  # For switches
    units::String
    description::String
end

"""
    introspect_actor_params(act::ParametersActors, actor_sym::Symbol)::Vector{ParamInfo}

Extract parameter information from an actor's parameters for UI rendering.

# Arguments
- `act::ParametersActors` - The full actor parameters struct
- `actor_sym::Symbol` - The actor name (e.g., :ActorFluxMatcher)

# Returns
- `Vector{ParamInfo}` - List of parameter information structs
"""
function introspect_actor_params(act::ParametersActors, actor_sym::Symbol)::Vector{ParamInfo}
    # Get actor parameter struct
    actor_params = getfield(act, actor_sym)

    params = ParamInfo[]

    # Iterate through fields
    for field_name in fieldnames(typeof(actor_params))
        # Skip internal fields (start with underscore)
        if startswith(string(field_name), "_")
            continue
        end

        field_val = getfield(actor_params, field_name)

        if field_val isa Entry
            # Extract Entry parameter info
            push!(params, ParamInfo(
                field_name,
                :entry,
                string(typeof(field_val).parameters[1]),
                field_val.value,
                field_val.default,
                [],  # No options for Entry
                field_val.units,
                field_val.description
            ))
        elseif field_val isa Switch
            # Extract Switch parameter info
            push!(params, ParamInfo(
                field_name,
                :switch,
                string(typeof(field_val).parameters[1]),
                field_val.value,
                field_val.default,
                collect(keys(field_val.options)),  # Available options
                field_val.units,
                field_val.description
            ))
        end
    end

    return params
end

"""
    update_param!(act::ParametersActors, actor_sym::Symbol, param_name::Symbol, new_value)

Update a parameter value in an actor's parameters.

# Arguments
- `act::ParametersActors` - The full actor parameters struct
- `actor_sym::Symbol` - The actor name
- `param_name::Symbol` - The parameter name
- `new_value` - The new value (as string, will be parsed)
"""
function update_param!(act::ParametersActors, actor_sym::Symbol, param_name::Symbol, new_value)
    actor_params = getfield(act, actor_sym)
    param = getfield(actor_params, param_name)

    if param isa Entry
        param.value = parse_param_value(new_value, typeof(param).parameters[1])
    elseif param isa Switch
        param.value = parse_param_value(new_value, typeof(param).parameters[1])
    else
        throw(ArgumentError("Parameter $param_name is not an Entry or Switch"))
    end
end

"""
    parse_param_value(str_value, ::Type{T}) where T

Parse a string value into the appropriate type.

# Arguments
- `str_value` - String representation of the value
- `T` - Target type

# Returns
- Parsed value of type T
"""
function parse_param_value(str_value, ::Type{T}) where {T}
    # Handle Symbol
    if T <: Symbol
        return Symbol(str_value)
    end

    # Handle Bool
    if T <: Bool
        lower_str = lowercase(string(str_value))
        if lower_str in ("true", "1", "yes", "y")
            return true
        elseif lower_str in ("false", "0", "no", "n")
            return false
        else
            throw(ArgumentError("Cannot parse '$str_value' as Bool"))
        end
    end

    # Handle numeric types
    if T <: Int
        return parse(Int, str_value)
    elseif T <: Float64
        return parse(Float64, str_value)
    end

    # Handle vectors (including ranges)
    if T <: AbstractVector
        # Try to evaluate as Julia expression (e.g., "0.25:0.1:0.85")
        try
            return eval(Meta.parse(str_value))
        catch e
            throw(ArgumentError("Cannot parse '$str_value' as $T: $e"))
        end
    end

    # Default: return as string
    return string(str_value)
end
