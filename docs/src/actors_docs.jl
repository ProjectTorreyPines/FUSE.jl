using InteractiveUtils: subtypes

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
  ActorEquilibrium  <--  generic
  ├─ ActorTEQUILA   <--  specific
  ├─ ActorCHEASE    <--  specific
  └─ ActorFRESCO    <--  specific
  ```
* `act.[ActorGeneric].model` selects specific actor being used
* All specific actors will expect data and fill the same enties in `dd`
  * IMAS.jl expressions are key to make this work seamlessly
* Where possible actors should make use of generic actors and not hardcode use of specific actors

```@contents
    Pages = ["actors.md"]
    Depth = 3
```
"""]

function concrete_subtypes(T::Type)
    if isabstracttype(T)
        sub = subtypes(T)
        return vcat(map(concrete_subtypes, sub)...)
    else
        return [T]
    end
end

list_directories(path::String) = (item for item in readdir(path) if isdir(joinpath(path, item)))

single_actors = concrete_subtypes(FUSE.SingleAbstractActor)
compound_actors = concrete_subtypes(FUSE.CompoundAbstractActor)

for actor_dir in list_directories(joinpath(FUSE.__FUSE__, "src", "actors"))

    first_time_actor_dir = true
    index_txt_header = 0
    n_actors = 0

    for Actor in [compound_actors; single_actors]
        folder = FUSE.group_name(Actor)
        name = FUSE.name(Actor; remove_Actor=false)
        nname = replace("$name", "Actor" => "")
        basename = replace(nname, "_" => " ")

        if folder == actor_dir
            if first_time_actor_dir
                push!(txt, "## $(uppercasefirst(replace(actor_dir,"_"=>" ")))")
                index_txt_header = length(txt)
                first_time_actor_dir = false
            end
            n_actors += 1
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
    if index_txt_header > 0
        txt[index_txt_header] = "$(txt[index_txt_header]) ($n_actors actors)"
    end
end

open(joinpath(@__DIR__, "actors.md"), "w") do io
    return write(io, join(txt, "\n"))
end