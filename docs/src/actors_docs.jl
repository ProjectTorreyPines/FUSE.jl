txt = ["""
# Physics and Engineering Actors

Physics and engineering **actors** are the fundamental building blocks of FUSE simulations:
  * Actors operate exclusively on `IMAS.dd` data
  * Actors functionality is controlled via `act` parameters
  * Actors can be combined into other actors
```@contents
    Pages = ["actors.md"]
    Depth = 3
```
"""]

for actor_abstract_type in subtypes(FUSE.AbstractActor)
    push!(txt, """## $(replace(replace("$actor_abstract_type","AbstractActor"=>""),"FUSE." => "")) actors\n""")
    for name in sort(collect(names(FUSE; all=true, imported=false)))
        if startswith("$name", "Actor") && supertype(@eval(FUSE, $name)) == actor_abstract_type
            nname = replace("$name", "Actor" => "")
            basename = replace(nname, "_" => " ")
            push!(txt,
                """### $basename

                ```@docs
                FUSE.$name(dd::IMAS.dd, act::FUSE.ParametersActor; kw...)
                ```

                ```@eval
                import Markdown, FUSE
                if !isempty(keys(FUSE.ParametersActor(:$name)))
                    return Markdown.parse("Valid `kw...` arguments are based on `FUSE.ParametersActor(:$name)`:")
                end
                ```

                ```@example
                import FUSE # hide
                act = FUSE.ParametersActor(:$name) # hide
                act._name=Symbol("act."*string(act._name)) # hide
                if !isempty(keys(act)) # hide
                    return act # hide
                end # hide
                ```
                """
            )
        end
    end
end

open("src/actors.md", "w") do io
    write(io, join(txt, "\n"))
end