txt = ["""
# ini Parameters

```@meta
CurrentModule = FUSE
```

```@example
import FUSE # hide
ini = FUSE.ParametersInit()
```

"""]
open("src/ini.md", "w") do io
    write(io, join(txt, "\n"))
end

open("src/ini_details.md", "w") do io
    parameters_details_md(io, FUSE.ParametersInit())
end