txt = ["""
# act Parameters

```@meta
CurrentModule = FUSE
```

```@example
import FUSE # hide
act = FUSE.ParametersActors() # hide
resize!(act.ActorSimpleEC.actuator, 1) # hide
resize!(act.ActorSimpleIC.actuator, 1) # hide
resize!(act.ActorSimpleLH.actuator, 1) # hide
resize!(act.ActorSimpleNB.actuator, 1) # hide
resize!(act.ActorSimplePL.actuator, 1) # hide
act # hide
```

"""]

open("$(@__DIR__)/act.md", "w") do io
    return write(io, join(txt, "\n"))
end

open("$(@__DIR__)/act_details.md", "w") do io
    act = FUSE.ParametersActors()
    resize!(act.ActorSimpleEC.actuator, 1) # hide
    resize!(act.ActorSimpleIC.actuator, 1) # hide
    resize!(act.ActorSimpleLH.actuator, 1) # hide
    resize!(act.ActorSimpleNB.actuator, 1) # hide
    resize!(act.ActorSimplePL.actuator, 1) # hide
    return parameters_details_md(io, act)
end
