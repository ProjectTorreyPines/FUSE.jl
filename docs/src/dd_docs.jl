txt = ["""
# IMAS data structure

FUSE data is organized into hierarchical Interface Data Structures (IDSs), according to the ITER IMAS ontology.
In addition to the usual IMAS IDSs (which we support on a need-by-need basis) FUSE also defines some of its own IDSs,
to hold data that does not (yet?) fit into IMAS. Notable examples are the `build`, `solid_mechanics`, `balance_of_plant`, and `costing` IDSs.

```@meta
CurrentModule = FUSE
```

## dd

`dd = IMAS.dd()` (which stands for "data dictionary") is the root of the FUSE data structure

```@example
using IMAS # hide
IMAS.dd # hide
```
"""]

for name in sort(collect(fieldnames(IMAS.dd)))
    if !startswith("$name", "_")
        basename = replace("$name", "_" => " ")
        push!(txt,
            """## $basename
            ```@example
            using IMASDD # hide
            IMASDD.$name # hide
            ```\n""")
    end
end
open("src/dd.md", "w") do io
    write(io, join(txt, "\n"))
end