"""
    insert_subtype_members T

Expands the subtypes into struct members
"""
macro insert_subtype_members(T)
    expr = Expr(:block)
    for s in subtypes(__module__.eval(T))
        full_parameters_actor_name = split(String(Symbol(s)), '.')[end]
        if !startswith(full_parameters_actor_name, '_')
            mem = Symbol(split(full_parameters_actor_name, "__")[end])
            # Splice the resolved subtype object directly into the AST. Using
            # an escaped symbol would defer the lookup to the caller module's
            # scope, requiring the caller to import every subtype by name —
            # which defeats the point of auto-discovering them via `subtypes`.
            push!(expr.args, Expr(:(::), esc(mem), Expr(:curly, s, esc(:T))))
        end
    end
    return expr
end

"""
    @insert_constructor_members T

Returns a tuple of all the evaluated cosntructors for the
subtypes we can then splat into the larger constructor
"""
macro insert_constructor_members(T)
    expr = Expr(:tuple)
    for s in subtypes(__module__.eval(T))
        full_actor_name = split(String(Symbol(s)), '.')[end]
        if !startswith(full_actor_name, '_')
            # Same hygiene rationale as @insert_subtype_members: splice the
            # resolved subtype object so the caller module does not have to
            # import each subtype's name.
            push!(expr.args, Expr(:call, Expr(:curly, s, esc(:T))))
        end
    end
    return expr
end

mutable struct ParametersActors{T<:Real} <: ParametersAllActors{T}
    _parent::WeakRef
    _name::Symbol
    @insert_subtype_members ParametersActor
end

function ParametersActors{T}() where {T<:Real}
    act = ParametersActors{T}(
        WeakRef(nothing),
        :act,
        (@insert_constructor_members ParametersActor)...
    )
    setup_parameters!(act)
    return act
end

function ParametersActors()
    return ParametersActors{Float64}()
end

############
# act_common_parameters
############
"""
    act_common_parameters(; kw...)

Returns commonly used act parameters as a switch or entry, example: act_common_parameters(do_plot=true)
"""
function act_common_parameters(; kw...)
    @assert length(kw) == 1 "act_common_parameters only takes one argument"
    name = first(keys(kw))
    default = first(values(kw))
    if name == :do_plot
        return Entry{Bool}("-", "Store the output dds of the workflow run"; default)
    elseif name == :verbose
        return Entry{Bool}("-", "Verbose"; default)
    else
        error("There is no act_common_parameter for name = $name")
    end
end

###############
# save / load #
###############

# Backward compatibility: rename old parameter keys to new names
const ACT_ACTOR_MIGRATIONS = Dict(
    "ActorHCD" => "ActorSources",
    "ActorSawteeth" => "ActorSawteethSource"
)

function _migrate_act_string(str::String)
    for (old, new) in ACT_ACTOR_MIGRATIONS
        str = replace(str, "\":$old\"" => "\":$new\"")
    end
    return str
end

function _assert_migrations_valid(act::ParametersAllActors)
    schema_keys = Set(string(k) for k in keys(act))
    for old_name in keys(ACT_ACTOR_MIGRATIONS)
        if old_name ∈ schema_keys
            error("Migration source '$old_name' conflicts with an existing actor — remove or update the entry in ACT_ACTOR_MIGRATIONS")
        end
    end
end

function _validate_act_keys(str::String, act::ParametersAllActors)
    known = Set("\":" * string(k) * "\"" for k in keys(act))
    unknown = String[]
    for m in eachmatch(r"\":(Actor\w+)\"", str)
        key = "\":$(m.captures[1])\""
        if key ∉ known
            push!(unknown, m.captures[1])
        end
    end
    unique!(unknown)
    if !isempty(unknown)
        error("Unknown actor(s) in act file: $(join(sort!(unknown), ", "))")
    end
end
"""
    act2json(act::ParametersAllActors, filename::AbstractString; kw...)

Save the ACT act parameters to a JSON file with given `filename`

`kw` arguments are passed to the JSON.print function
"""
function act2json(act::ParametersAllActors, filename::AbstractString; kw...)
    return SimulationParameters.par2json(act, filename; kw...)
end

"""
    json2act(filename::AbstractString, act::ParametersAllActors=ParametersActors())

Load the ACT act parameters from a JSON file with given `filename`
"""
function json2act(filename::AbstractString, act::ParametersAllActors=ParametersActors())
    _assert_migrations_valid(act)
    str = _migrate_act_string(read(filename, String))
    _validate_act_keys(str, act)
    return SimulationParameters.jstr2par(str, act)
end

"""
    act2yaml(act::ParametersAllActors, filename::AbstractString; kw...)

Save the ACT parameters to a YAML file with given `filename`

`kw` arguments are passed to the YAML.print function
"""
function act2yaml(act::ParametersAllActors, filename::AbstractString; kw...)
    return SimulationParameters.par2yaml(act, filename; kw...)
end

"""
    yaml2act(filename::AbstractString, act::ParametersAllActors=ParametersActors())

Load the ACT act parameters from a YAML file with given `filename`
"""
function yaml2act(filename::AbstractString, act::ParametersAllActors=ParametersActors())
    _assert_migrations_valid(act)
    str = _migrate_act_string(read(filename, String))
    _validate_act_keys(str, act)
    return SimulationParameters.ystr2par(str, act)
end

"""
    act2dict(act::ParametersAllInits; kw...)

Convert the ACT parameters to a dictionary form
"""
function act2dict(act::ParametersAllActors; kw...)
    return SimulationParameters.par2dict(act; kw...)
end

"""
    dict2act(dict::AbstractDict, act::ParametersAllInits=ParametersActors())

Convert dict to ACT parameters
"""
function dict2act(dict::AbstractDict, act::ParametersAllActors=ParametersActors())
    return SimulationParameters.dict2par!(dict, act)
end