import YAML
import OrderedCollections

function yaml_workflow(steps::Union{AbstractDict,AbstractVector}, dd::IMAS.dd, ini::ParametersInits, act::ParametersActors)
    actors = AbstractActor[]
    for step in steps
        if typeof(step) <: AbstractDict
            step_name, kwargs = first(step)
            step_name = Symbol(step_name)
            kwargs = SimulationParameters.replace_colon_strings_to_symbols(kwargs)
        else
            step_name = Symbol(step)
            kwargs = Dict()
        end
        if startswith(string(step_name), "init")
            getfield(FUSE, step_name)(dd, ini, act)
        elseif startswith(string(step_name), "Actor")
            actor = getfield(FUSE, step_name)(dd, act; kwargs...)
            push!(actors, actor)
        else
            @error("Invalid step `$(step_name)`")
        end
    end
    return actors
end

function yaml_workflow(flow::String, dd::IMAS.dd, ini::ParametersInits, act::ParametersActors)
    steps = YAML.load(flow; dicttype=OrderedCollections.OrderedDict)
    return yaml_workflow(steps, dd, ini, act)
end
