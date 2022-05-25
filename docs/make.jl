using Documenter
import FUSE
import IMAS
import IMASDD
import AbstractTrees

function html_link_repr(par::FUSE.Parameter)
    return "©" * join(FUSE.path(par), ".") * "©©" * string(par._name) * "©"
end
function AbstractTrees.printnode(io::IO, par::FUSE.Parameter)
    return printstyled(io, html_link_repr(par))
end

function parameters_details_md(io, pars)
    for leaf in AbstractTrees.Leaves(pars)
        if typeof(leaf) <: FUSE.Parameters
            continue
        end
        if ismissing(leaf.default)
            default = ""
            note = "note"
        else
            default = "* **Default:** `$(leaf.default)`"
            note = "tip"
        end
        if typeof(leaf) <: FUSE.Switch
            options = "* **Options:** " * join(["`$(opt.first)`" for opt in leaf.options], ", ")
        else
            options = ""
        end
        txt = """

        ------------

        ```@raw html
        <div id='$(join(FUSE.path(leaf),"."))'>
        ```
        ## `$(join(FUSE.path(leaf),"."))`

        !!! $note "$(leaf.description)"
            * **Units:** `$(isempty(leaf.units) ? "-" : leaf.units)`
            $options
            $default

        """
        write(io, txt)
    end
end

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
makedocs(;
    modules=[FUSE, IMAS, IMASDD],
    sitename="FUSE",
    format=Documenter.HTML(; prettyurls=false, sidebar_sitename=false, assets=["assets/favicon.ico"]),
    pages=[
        "index.md",
        "dd Data Structure" => "dd.md",
        "ini Parameters" => "ini.md",
        "act Parameters" => "act.md",
        "Actors" => "actors.md",
        "Initialization" => "inits.md",
        "Use Cases" => "cases.md",
        "Utilities" => "utils.md",
        "Installation" => "install.md",
    ],
)

# convert "©(.*)©©(.*)©" patterns to hyperlinks
for (file, parfile) in [("act", "act"), ("ini", "ini"), ("actors", "act")]
    open("build/$file.html", "r") do io
        txt = read(io, String)
        txt = split(txt, "\n")
        for (k, line) in enumerate(txt)
            txt[k] = replace(replace(line, r"©(.*)©©(.*)©" => s"<a href='©_details.html#\1'>\2</a>"), "©" => parfile)
        end
        open("build/$file.html", "w") do io
            write(io, join(txt, "\n"))
        end
    end
end

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
