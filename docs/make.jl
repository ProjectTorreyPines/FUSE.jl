using Documenter
import FUSE, IMAS, IMASDD

# ================ #
# generate DD page #
# ================ #
include("src/dd_docs.jl")

# ==================== #
# generate Actors page #
# ==================== #
include("src/actors_docs.jl")

# =================== #
# generate inits page #
# =================== #
include("src/inits_docs.jl")

# ================= #
# generate ini page #
# ================= #
include("src/ini_docs.jl")

# ================= #
# generate act page #
# ================= #
include("src/act_docs.jl")

# =================== #
# generate cases page #
# =================== #
include("src/cases_docs.jl")

# ============== #
# build the docs #
# ============== #
makedocs(
    modules=[FUSE,IMAS,IMASDD],
    sitename="FUSE",
    format=Documenter.HTML(prettyurls=false, sidebar_sitename=false, assets=["assets/favicon.ico"]),
    pages = [
        "index.md",
        "dd Data Structure" => "dd.md",
        "ini Parameters" => "ini.md",
        "act Parameters" => "act.md",
        "Actors" => "actors.md",
        "Initialization" => "inits.md",
        "Use Cases" => "cases.md",
        "Utilities" => "utils.md",
        "Installation" => "install.md"
        ]
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
