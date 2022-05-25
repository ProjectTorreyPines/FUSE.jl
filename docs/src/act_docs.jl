txt = ["""
# Actor Parameters

```@meta
CurrentModule = FUSE
```

```@example
import FUSE # hide
act = FUSE.ParametersActor()
```

"""]

open("src/act.md", "w") do io
    write(io, join(txt, "\n"))
end

open("src/act_details.md", "w") do io
    parameters_details_md(io, FUSE.ParametersActor())
end