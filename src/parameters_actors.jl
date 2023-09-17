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


mutable struct ParametersActors{T} <: ParametersAllActors where {T<:Real}
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

"""
    act2json(act::ParametersAllActors, filename::AbstractString; kw...)

Save the FUSE parameters to a JSON file with give `filename`
`kw` arguments are passed to the JSON.print function
"""
function act2json(act::ParametersAllActors, filename::AbstractString; kw...)
    return SimulationParameters.par2json(act, filename; kw...)
end

function json2act(filename::AbstractString)
    return SimulationParameters.json2par(filename, ParametersActors())
end