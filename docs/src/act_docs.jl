txt = ["""
# act Parameters

```@meta
CurrentModule = FUSE
```

```@example
import FUSE # hide
act = FUSE.ParametersAllActors()
```

"""]

open("$(@__DIR__)/act.md", "w") do io
    write(io, join(txt, "\n"))
end

open("$(@__DIR__)/act_details.md", "w") do io
    parameters_details_md(io, FUSE.ParametersAllActors())
end