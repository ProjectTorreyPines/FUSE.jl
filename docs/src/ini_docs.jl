txt = ["""
# ini Parameters

```@meta
CurrentModule = FUSE
```

```@example
import FUSE # hide
ini = FUSE.ParametersAllInits()
```

"""]
open("$(@__DIR__)/ini.md", "w") do io
    write(io, join(txt, "\n"))
end

open("$(@__DIR__)/ini_details.md", "w") do io
    parameters_details_md(io, FUSE.ParametersAllInits())
end