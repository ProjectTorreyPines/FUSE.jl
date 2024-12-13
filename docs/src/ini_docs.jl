txt = ["""
# ini Parameters

```@meta
CurrentModule = FUSE
```

```@example
import FUSE # hide
ini = FUSE.ParametersInits() # hide
resize!(ini.build.layers, 1) # hide
resize!(ini.nb_unit, 1) # hide
resize!(ini.ec_launcher, 1) # hide
resize!(ini.pellet_launcher, 1) # hide
resize!(ini.ic_antenna, 1) # hide
resize!(ini.lh_antenna, 1) # hide
ini # hide
```

"""]
open("$(@__DIR__)/ini.md", "w") do io
    return write(io, join(txt, "\n"))
end

open("$(@__DIR__)/ini_details.md", "w") do io
    ini = FUSE.ParametersInits()
    resize!(ini.build.layers, 1)
    resize!(ini.nb_unit, 1)
    resize!(ini.ec_launcher, 1)
    resize!(ini.pellet_launcher, 1)
    resize!(ini.ic_antenna, 1)
    resize!(ini.lh_antenna, 1)
    return parameters_details_md(io, ini)
end