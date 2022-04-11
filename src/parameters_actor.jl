using InteractiveUtils: subtypes
"""
    ActorParameters()

Generates actor parameters 
"""
function ActorParameters()
    act = ActorParameters(Symbol[], Dict{Symbol,Union{Parameter,ActorParameters}}())
    for par in subtypes(FUSE.AbstractActor)
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

