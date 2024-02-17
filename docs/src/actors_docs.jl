txt = ["""
# Physics and Engineering Actors

Physics and engineering **actors** are the fundamental building blocks of FUSE simulations:
* Actors operate exclusively on `IMAS.dd` data
* Actors functionality is controlled via `act` parameters
* Actors can be combined into other actors

Fidelity hierarchy is enabled by concept of *generic* Vs *specific* actors
* Generic actors define physics/component 
* Specific actors implement a specific model for that physics/component
* For example:
  ```
  EquilibriumActor  <--  generic
  ├─ SolovevActor   <--  specific
  └─ CHEASEActor    <--  specific
  ```
* `act.[GenericActor].model` selects specific actor being used
* All specific actors will expect data and fill the same enties in `dd`
  * IMAS.jl expressions are key to make this work seamlessly
* Where possible workflows should make use of generic actors and not hardcode use of specific actors

```@contents
    Pages = ["actors.md"]
    Depth = 3
```
"""]

list_directories(path::String) = [item for item in readdir(path) if isdir(joinpath(path, item))]
act = FUSE.ParametersActors()
for actor_dir in list_directories(joinpath(FUSE.__FUSE__, "src", "actors"))
    push!(txt, """## $actor_dir actors""")
    for name in sort!(collect(names(FUSE; all=true, imported=false)))
        if !startswith("$name", "Actor") #|| name ∈ skip_list
            continue
        end
        folder = ""
        if supertype(@eval(FUSE, $name)) == FUSE.SingleAbstractActor
            which_output = string(@which getfield(FUSE, name)(IMAS.dd(), getfield(act, name)))
            folder = split(split(which_output, "@")[end], "/")[end-1]
        elseif supertype(@eval(FUSE, $name)) == FUSE.CompoundAbstractActor
            which_output = string(@which getfield(FUSE, name)(IMAS.dd(), getfield(act, name), act))
            folder = split(split(which_output, "@")[end], "/")[end-1]
        end
        if !isempty(folder) && folder == actor_dir
            nname = replace("$name", "Actor" => "")
            basename = replace(nname, "_" => " ")
            push!(txt,
                """### $basename

                ```@docs
                FUSE.$name(dd::IMAS.dd, act::FUSE.ParametersAllActors; kw...)
                ```

                ```@example
                import FUSE # hide
                act = FUSE.ParametersActors() # hide
                getfield(FUSE.ParametersActors(), :$name) # hide
                ```
                """
            )
        end
    end
end

open("$(@__DIR__)/actors.md", "w") do io
    return write(io, join(txt, "\n"))
end