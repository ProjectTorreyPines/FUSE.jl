"""
    insert_subtype_members T

Expands the subtypes into struct members
"""
macro insert_subtype_members(T)
    expr = Expr(:block)
    for s in subtypes(eval(T))
        mem = Symbol(replace(String(Symbol(s)), "FUSE.FUSEparameters__" => ""))
        constraint = Symbol(replace(String(Symbol(s)), "FUSE." => ""))
        push!(expr.args, Expr(:(::), esc(mem), Expr(:curly, esc(constraint), esc(:T))))
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
    for s in subtypes(eval(T))
        mem = Symbol(replace(String(Symbol(s)), "FUSE." => ""))
        push!(expr.args, Expr(:call, Expr(:curly, esc(mem), esc(:T))))
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
        return  Entry{Bool}("-", "Verbose"; default)
    else
        error("There is no act_common_parameter for name = $name")
    end
end

###############
# save / load #
###############

"""
    act2json(act::ParametersAllActors, filename::AbstractString; kw...)

Save the FUSE act parameters to a JSON file with given `filename`

`kw` arguments are passed to the JSON.print function
"""
function act2json(act::ParametersAllActors, filename::AbstractString; kw...)
    return SimulationParameters.par2json(act, filename; kw...)
end

"""
    json2act(filename::AbstractString, act::ParametersAllActors=ParametersActors())

Load the FUSE act parameters from a JSON file with given `filename`
"""
function json2act(filename::AbstractString, act::ParametersAllActors=ParametersActors())
    return SimulationParameters.json2par(filename, act)
end

"""
    act2yaml(act::ParametersAllActors, filename::AbstractString; kw...)

Save the FUSE parameters to a YAML file with given `filename`

`kw` arguments are passed to the YAML.print function
"""
function act2yaml(act::ParametersAllActors, filename::AbstractString; kw...)
    return SimulationParameters.par2yaml(act, filename; kw...)
end

"""
    yaml2act(filename::AbstractString, act::ParametersAllActors=ParametersActors())

Load the FUSE act parameters from a YAML file with given `filename`
"""
function yaml2act(filename::AbstractString, act::ParametersAllActors=ParametersActors())
    return SimulationParameters.yaml2par(filename, act)
end
