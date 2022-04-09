"""
    ActorParameters()

Generates actor parameters 
"""
function ActorParameters()
    par = ActorParameters(Symbol[], Dict{Symbol,Union{Parameter,ActorParameters}}())
    for item in [:SolovevActor]
        setproperty!(par, item, ActorParameters(item))
    end
    return par
end

function ActorParameters(::Type{Val{:SolovevActor}})
    par = ActorParameters(nothing)
    par.ngrid = Entry(Integer, "", "ngrid"; default=129)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end
