using Pkg
Pkg.activate("..")

using Revise

using Documenter
import FUSE
import IMAS
import IMASDD
import AbstractTrees
import ProgressMeter

function html_link_repr(par::FUSE.AbstractParameter)
    return "©" * join(FUSE.path(par), ".") * "©©" * string(par._name) * "©"
end

function AbstractTrees.printnode(io::IO, par::FUSE.AbstractParameter)
    return printstyled(io, html_link_repr(par))
end

function AbstractTrees.printnode(io::IO, @nospecialize(ids::Type{T}); kwargs...) where {T<:IMAS.IDS}
    txt = split(split("$ids", "___")[end], "__")[end]
    return printstyled(io, txt; bold=true)
end

function AbstractTrees.printnode(io::IO, leaf::IMAS.IMASleafRepr; kwargs...)
    if startswith(leaf.location, "dd.")
        printstyled(io, "©$(replace(leaf.location[4:end],"_"=>"-"))©©©$(leaf.key)©")
    else
        printstyled(io, "©$(leaf.location)©©$(leaf.key)©")
    end
end

function parameters_details_md(io::IO, pars::FUSE.AbstractParameters)
    for leaf in AbstractTrees.Leaves(pars)
        if typeof(leaf) <: FUSE.AbstractParameters
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
            options = "* **Options:** " * join(["`$(opt.first)`" for opt in leaf.options], ", ") * "\n    "
        else
            options = ""
        end
        txt = """

        ------------

        ```@raw html
        <div id='$(join(FUSE.path(leaf),"."))'></div>
        ```
        !!! $note "$(join(FUSE.path(leaf),"."))"
            $(leaf.description)
            * **Units:** `$(isempty(leaf.units) ? "-" : leaf.units)`
            $(options)$(default)

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

# ====================== #
# generate examples page #
# ====================== #
include("src/examples.jl")

# ============== #
# build the docs #
# ============== #
makedocs(;
    modules=[FUSE, IMAS, IMASDD],
    sitename="FUSE",
    format=Documenter.HTML(; prettyurls=false, sidebar_sitename=false, assets=["assets/favicon.ico"]),
    pages=[
        "Concepts" => "index.md",
        "Data Structure" => "dd.md",
        "Actors" => "actors.md",
        "Parameters" => ["ini Parameters" => "ini.md", "act Parameters" => "act.md", "Use Cases" => "cases.md", "Initialization" => "inits.md"],
        "Examples" => "examples.md",
        "Development" => "develop.md",
        "Install" => "install.md",
        "Others" => ["GASC" => "gasc.md", "Utilities" => "utils.md", "HPC" => "parallel.md"],
    ],
)

# convert "©(.*)©©(.*)©" patterns to hyperlinks
@info "Converting links"
for (file, parfile) in [("act", "act"), ("ini", "ini"), ("actors", "act"), ("dd", "dd")]
    local txt = open("build/$file.html", "r") do io
        read(io, String)
    end
    txt = split(txt, "\n")
    for (k, line) in enumerate(txt)
        txt[k] = replace(txt[k], r"©(.*)©©©(.*)©" => s"<a href='©.html#\1'>\2</a>")
        txt[k] = replace(txt[k], r"©(.*)©©(.*)©" => s"<a href='©_details.html#\1'>\2</a>")
        txt[k] = replace(txt[k], "©" => parfile)
    end
    open("build/$file.html", "w") do io
        write(io, join(txt, "\n"))
    end
end

# distinguish between input/output cells
@info "Styling examples"
for css in ["light", "dark"]
    open("build/assets/themes/documenter-$css.css", "a") do io
        write(
            io,
            """\n
.nohighlight {
background-color: transparent !important;
}""",
        )
    end
end
files_to_convert = readdir("build")[findall(x -> startswith(x, "example_") && endswith(x, ".html"), readdir("build"))]
for file in files_to_convert
    local txt = open("build/$file", "r") do io
        read(io, String)
    end
    txt = split(txt, "\n")
    for (k, line) in enumerate(txt)
        txt[k] = replace(line, "<pre><code class=\"nohighlight" => "<pre class=\"nohighlight\"><code class=\"nohighlight")
    end
    open("build/$file", "w") do io
        write(io, join(txt, "\n"))
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
