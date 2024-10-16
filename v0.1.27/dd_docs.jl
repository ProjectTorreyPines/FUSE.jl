txt = String["""

```@meta
CurrentModule = IMAS
```

# IMAS data structure

FUSE data is organized into hierarchical Interface Data Structures (IDSs), according to the ITER IMAS ontology.
In addition to the usual IMAS IDSs (which we include on a need-by-need basis) FUSE also defines some of its own IDSs,
to hold data that does not (yet?) fit into IMAS. Notable examples are the `build`, `solid_mechanics`, `balance_of_plant`, and `costing` IDSs.

`dd = IMAS.dd()` (which stands for "data dictionary") is the root of the FUSE data structure

"""]

for name in sort!(collect(fieldnames(IMAS.dd)))
    if name == :global_time || name ∈ IMAS.private_fields
        continue
    else
        basename = replace("$name", "_" => " ")
        push!(
            txt,
            """## $basename
            ```@example
            using IMASdd # hide
            IMASdd.$name{Float64} # hide
            ```\n"""
        )
    end
end
open("$(@__DIR__)/dd.md", "w") do io
    return write(io, join(txt, "\n"))
end

function dd_details_md(io, ids)
    ProgressMeter.@showprogress "$ids" for leaf in collect(AbstractTrees.Leaves(ids))
        name = "$(leaf.location)"
        nfo = IMAS.info(name)
        documentation = nfo.documentation
        if nfo.documentation != ""
            documentation = nfo.documentation
        else
            documentation = "N/A"
        end
        if nfo.units != "-"
            units = "* **Units:** `$(nfo.units)`\n    "
        else
            units = ""
        end
        data_type = "* **Data Type:** `$(nfo.data_type)`\n    "
        if !isempty(nfo.coordinates)
            coordinates = "* **Coordinates:** `$(String[k for k in nfo.coordinates])`\n    "
        else
            coordinates = ""
        end
        txt = """

        ------------

        ```@raw html
        <div id='$name'></div>
        ```
        !!! note "$name"
            $documentation
            $(units)$(data_type)$(coordinates)
        """
        write(io, txt)
    end
end

open("$(@__DIR__)/dd_details.md", "w") do io
    for field in fieldnames(IMAS.dd)
        if field == :global_time || field ∈ IMAS.private_fields
            continue
        end
        dd_details_md(io, getfield(IMAS, field){Float64})
    end
end
