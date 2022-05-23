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
        push!(txt,"""## $(replace("$name","_"=>" "))
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
        basename = replace("$name", "Actor"=>"")
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
    if startswith("$name", "init")
        basename = replace("$name", "init_"=>"")
        push!(
            txt,
            """## $basename

```@docs
FUSE.$name(::IMAS.dd, ::FUSE.ParametersInit, ::FUSE.ParametersActor)
```

```@eval
import Markdown, FUSE
return Markdown.parse(FUSE.doc(FUSE.ParametersInit(:$basename)))
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
