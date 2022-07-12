file_dir = dirname(abspath(@__FILE__))
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

open("$file_dir/act.md", "w") do io
    write(io, join(txt, "\n"))
end

open("$file_dir/act_details.md", "w") do io
    parameters_details_md(io, FUSE.ParametersAllActors())
end