txt = ["""
# act Parameters

```@meta
CurrentModule = FUSE
```

```@example
import FUSE # hide
act = FUSE.ParametersActors()
```

"""]

open("$(@__DIR__)/act.md", "w") do io
    return write(io, join(txt, "\n"))
end

open("$(@__DIR__)/act_details.md", "w") do io
    return parameters_details_md(io, FUSE.ParametersActors())
end