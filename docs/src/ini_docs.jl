file_dir = dirname(abspath(@__FILE__))
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
open("$file_dir/ini.md", "w") do io
    write(io, join(txt, "\n"))
end

open("$file_dir/ini_details.md", "w") do io
    parameters_details_md(io, FUSE.ParametersAllInits())
end