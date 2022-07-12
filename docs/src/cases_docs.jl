file_dir = dirname(abspath(@__FILE__))
txt = ["""# Use cases

```@meta
CurrentModule = FUSE
```

FUSE comes with a set of pre-cookes used cases.
The `case_parameters(:use_case, ...)` method returns the `ini` and `act` parameters for that specific `use_case`.
These `ini` and `act` can then be further customized before running a FUSE simulation.

!!! tip "Tip!"
    Click on the `Source` button of each use case to see how each is setup

"""]
for method in methods(FUSE.case_parameters)
    name = try
        method.sig.types[2].parameters[1].parameters[1]
    catch
        continue
    end
    push!(txt,
        """## $name

        ```@docs
        case_parameters(::Type{Val{:$name}})
        ```
        """
    )
end
open("$file_dir/cases.md", "w") do io
    write(io, join(txt, "\n"))
end