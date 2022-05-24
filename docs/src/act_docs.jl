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