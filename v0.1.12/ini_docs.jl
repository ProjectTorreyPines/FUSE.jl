txt = ["""
# ini Parameters

```@meta
CurrentModule = FUSE
```

```@example
import FUSE # hide
ini = FUSE.ParametersInits(; n_ic=1, n_nb=1, n_ec=1, n_lh=1, n_pl=1, n_layers=1)
```

"""]
open("$(@__DIR__)/ini.md", "w") do io
    return write(io, join(txt, "\n"))
end

open("$(@__DIR__)/ini_details.md", "w") do io
    return parameters_details_md(io, FUSE.ParametersInits(; n_ic=1, n_nb=1, n_ec=1, n_lh=1, n_pl=1, n_layers=1))
end