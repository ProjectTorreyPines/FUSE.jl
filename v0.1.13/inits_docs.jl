txt = ["""
# Initialization

```@meta
CurrentModule = IMAS
```
Recall that FUSE actors operate exclusively on `IMAS.dd` data.
As such, to run any actor, one has to first **initialize** `IMAS.dd` with some data.
This can be done in a number of ways:

  1. Manually/Interactivelly (eg. in the REPL or a Jupyter sesion)
  1. Starting from 0D [`ini`](ini.md) and [`act`](ini.md) parameters (same spirit of OMFIT's PRO_create module)
  1. By reading in an existing OMAS JSON data structure with [`json2imas`](@ref)
  1. Starting from [GA Systems Code](@ref) output, then to `ini`, and finally to `dd`

The following `init...()` routines initialize `dd` from 0D parameters (method #2)

```@meta
CurrentModule = FUSE
```

## High-level Initialization

```@docs
FUSE.init(::IMAS.dd, ::FUSE.ParametersAllInits, ::FUSE.ParametersAllActors)
```

## Use-cases initialization

```@docs
FUSE.init(case::Symbol; do_plot:Bool=false, kw...)
```

-----------

## Low-level initialization routines

Below are the initialization functions specific to IDSs in the `dd` data structure.
These can be called for a fine control on what IDSs are initialized and how.

"""]
for name in sort!(collect(names(FUSE; all=true, imported=false)))
    if startswith("$name", "init_")
        nname = replace("$name", "init_" => "")
        basename = replace(nname, "_" => " ")
        push!(txt,
            """### $basename

            ```@docs
            FUSE.$name(::IMAS.dd, ::FUSE.ParametersAllInits, ::FUSE.ParametersAllActors)
            ```

            """
        )
    end
end
open("$(@__DIR__)/inits.md", "w") do io
    return write(io, join(txt, "\n"))
end
