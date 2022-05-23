using Documenter, FUSE

# ================ #
# generate DD page #
# ================ #
txt = ["""
# Data dictionary

## dd

```@example
using IMASDD # hide
IMASDD.dd # hide
```
"""]
for name in fieldnames(IMAS.dd)
    if !startswith("$name","_")
        basename = replace("$name","_"=>" ")
        push!(txt,"""## $basename
```@example
using IMASDD # hide
IMASDD.$name # hide
```\n""")
    end
end
open("src/dd.md", "w") do io
    write(io, join(txt, "\n"))
end

# ==================== #
# generate Actors page #
# ==================== #
txt = ["# Actors\n"]
for name in names(FUSE; all=true, imported=false)
    if startswith("$name", "Actor")
        nname=replace("$name", "Actor"=>"")
        basename = replace(nname,"_"=>" ")
        push!(
            txt,
            """## $basename

```@docs
FUSE.$name(::IMAS.dd, ::FUSE.ParametersActor)
```

```@eval
import Markdown, FUSE
return Markdown.parse(FUSE.doc(FUSE.ParametersActor(:$name)))
```
"""
        )
    end
end
open("src/actors.md", "w") do io
    write(io, join(txt, "\n"))
end

# ================= #
# generate init page #
# ================= #
txt = ["# Init\n"]
for name in names(FUSE; all=true, imported=false)
    if startswith("$name", "init_")
        nname = replace("$name", "init_"=>"")
        basename = replace(nname,"_"=>" ")
        push!(
            txt,
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

# ============== #
# build the docs #
# ============== #
makedocs(
    modules=[FUSE, IMAS],
    sitename="FUSE",
    format=Documenter.HTML(prettyurls=false)
)

# # =============== #
# # deploy the docs #
# # =============== #
# deploydocs(
#     target = "build",
#     repo = "github.com:ProjectTorreyPines/FUSE.jl.git",
#     forcepush = true
# )

# makedocs(
#     modules=[FUSE, IMAS],
#     sitename="FUSE",
#     format=Documenter.HTML(prettyurls=false)
# )
