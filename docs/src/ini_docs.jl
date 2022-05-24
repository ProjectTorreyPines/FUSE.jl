txt = ["""
# Init Parameters

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