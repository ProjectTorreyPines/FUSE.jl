txt = ["""
# Initialization

```@meta
CurrentModule = FUSE
```
Recall that FUSE actors operate exclusively on `IMAS.dd` data.
As such, to run any actor, one has to first **initialize** `IMAS.dd` with some data.
This can be done

  * Manually/Interactivelly (eg. in the REPL or a Jupyter sesion)
  * By reading in an existing OMAS JSON data structure
  * Starting from 0D `ini` parameters (same spirit of OMFIT's PRO_create module)
  * Starting from GASC outputs, then to `ini`, and finally to `dd`

## High level initialization

```@docs
FUSE.init(::IMAS.dd, ::FUSE.ParametersInit, ::FUSE.ParametersActor)
```

## Use-cases initialization

```@docs
FUSE.init(case::Symbol; do_plot:Bool=false, kw...)
```

"""]
for name in sort(collect(names(FUSE; all=true, imported=false)))
    if startswith("$name", "init_")
        nname = replace("$name", "init_" => "")
        basename = replace(nname, "_" => " ")
        push!(txt,
            """## $basename

            ```@docs
            FUSE.$name(::IMAS.dd, ::FUSE.ParametersInit, ::FUSE.ParametersActor)
            ```

            ```@eval
            import Markdown, FUSE
            return Markdown.parse(FUSE.doc(FUSE.ParametersInit(:$nname)))
            ```
            """
        )
    end
end
open("src/inits.md", "w") do io
    write(io, join(txt, "\n"))
end
