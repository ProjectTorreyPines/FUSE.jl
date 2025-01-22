function consistent_ini_act!(ini::ParametersAllInits, act::ParametersAllActors)
    if !isempty(ini.ec_launcher)
        if isempty(act.ActorSimpleEC.actuator)
            resize!(act.ActorSimpleEC, length(ini.ec_launcher))
        else
            @assert length(act.ActorSimpleEC.actuator) == length(ini.ec_launcher) "length(act.ActorSimpleEC.actuator) = $(length(act.ActorSimpleEC.actuator)) must be equal to length(ini.ec_launcher)=$(length(ini.ec_launcher))"
        end
    end
    if !isempty(ini.ic_antenna)
        if isempty(act.ActorSimpleIC.actuator)
            resize!(act.ActorSimpleIC, length(ini.ic_antenna))
        else
            @assert length(act.ActorSimpleIC.actuator) == length(ini.ic_antenna) "length(act.ActorSimpleIC.actuator) = $(length(act.ActorSimpleIC.actuator)) must be equal to length(ini.ic_antenna)=$(length(ini.ic_antenna))"
        end
    end
    if !isempty(ini.lh_antenna)
        if isempty(act.ActorSimpleLH.actuator)
            resize!(act.ActorSimpleLH, length(ini.lh_antenna))
        else
            @assert length(act.ActorSimpleLH.actuator) == length(ini.lh_antenna) "length(act.ActorSimpleLH.actuator) = $(length(act.ActorSimpleLH.actuator)) must be equal to length(ini.lh_antenna)=$(length(ini.lh_antenna))"
        end
    end
    if !isempty(ini.nb_unit)
        if isempty(act.ActorSimpleNB.actuator)
            resize!(act.ActorSimpleNB, length(ini.nb_unit))
        else
            @assert length(act.ActorSimpleNB.actuator) == length(ini.nb_unit) "length(act.ActorSimpleNB.actuator) = $(length(act.ActorSimpleNB.actuator)) must be equal to length(ini.nb_unit)=$(length(ini.nb_unit))"
        end
    end
    if !isempty(ini.pellet_launcher)
        if isempty(act.ActorSimplePellet.actuator)
            resize!(act.ActorSimplePellet, length(ini.pellet_launcher))
        else
            @assert length(act.ActorSimplePellet.actuator) == length(ini.pellet_launcher) "length(act.ActorSimplePellet.actuator) = $(length(act.ActorSimplePellet.actuator)) must be equal to length(ini.pellet_launcher)=$(length(ini.pellet_launcher))"
        end
    end
end

"""
    insert_subtype_members T

Expands the subtypes into struct members
"""
macro insert_subtype_members(T)
    expr = Expr(:block)
    for s in subtypes(eval(T))
        if !startswith(String(Symbol(s)), "FUSE._")
            mem = Symbol(replace(String(Symbol(s)), "FUSE.FUSEparameters__" => ""))
            constraint = Symbol(replace(String(Symbol(s)), "FUSE." => ""))
            push!(expr.args, Expr(:(::), esc(mem), Expr(:curly, esc(constraint), esc(:T))))
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
    for s in subtypes(eval(T))
        if !startswith(String(Symbol(s)), "FUSE._")
            mem = Symbol(replace(String(Symbol(s)), "FUSE." => ""))
            push!(expr.args, Expr(:call, Expr(:curly, esc(mem), esc(:T))))
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
    return SimulationParameters.json2par(filename, act)
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
    return SimulationParameters.yaml2par(filename, act)
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