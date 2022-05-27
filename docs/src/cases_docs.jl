txt = ["""# Use cases

```@meta
CurrentModule = FUSE
```

FUSE comes with a set of pre-cookes used cases. These use-cases are stored under the `FUSE/cases` folder.
The `case_parameters(:use_case, ...)` method returns the `ini` and `act` parameters for that specific `use_case`.

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
open("src/cases.md", "w") do io
    write(io, join(txt, "\n"))
end