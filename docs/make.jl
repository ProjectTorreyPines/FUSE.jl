using Pkg
Pkg.activate(@__DIR__)
using Documenter
import FUSE
import IMAS
import IMASdd
import SimulationParameters
import AbstractTrees
import ProgressMeter
import Dates
using Literate

include("notebook_to_jl.jl")

# Convert the Jupyter Notebook to Literate script
if isfile(joinpath(@__DIR__, "..", "examples", "tutorial.ipynb"))
    convert_notebook_to_litterate(joinpath(@__DIR__, "..", "examples", "tutorial.ipynb"), joinpath(@__DIR__, "src", "tutorial.jl"))
end
# Convert the Literate script to markdown
Literate.markdown(joinpath(@__DIR__, "src", "tutorial.jl"), joinpath(@__DIR__, "src"); documenter=true)
lines = join(readlines(joinpath(@__DIR__, "src", "tutorial.md")), "\n")
open(joinpath(@__DIR__, "src", "tutorial.md"), "w") do f
    return write(f, replace(lines,
        """
        ---

        *This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
        """ => ""))
end

function pretty_units(unit)
    unit = replace(unit, r"\^-3(?![0-9])" => "⁻³")
    unit = replace(unit, r"\^-2(?![0-9])" => "⁻²")
    unit = replace(unit, r"\^-1(?![0-9])" => "⁻¹")
    unit = replace(unit, r"\^3(?![0-9])" => "³")
    unit = replace(unit, r"\^2(?![0-9])" => "²")
    return unit
end

function html_link_repr(field::Symbol, location::AbstractString)
    return html_link_repr(string(field), location)
end

function html_link_repr(field::AbstractString, location::AbstractString)
    return "©" * location * "©©" * field * "©"
end

function AbstractTrees.printnode(io::IO, node_value::SimulationParameters.ParsNodeRepr)
    field = node_value.field
    par = node_value.value
    location = join(SimulationParameters.path(par), ".")
    if typeof(par) <: SimulationParameters.AbstractParameters
        printstyled(io, field; bold=true)
    elseif typeof(par) <: SimulationParameters.AbstractParametersVector
        printstyled(io, field; bold=true)
    elseif typeof(par) <: SimulationParameters.AbstractParameter
        if typeof(par.value) <: AbstractDict
            printstyled(io, "$field[:]"; bold=true)
        else
            units = par.units
            if units == "-"
                units = ""
            else
                units = "  [$(pretty_units(units))]"
            end
            printstyled(io, "$(html_link_repr(field, location))$units")
        end
    end
end

function AbstractTrees.printnode(io::IO, @nospecialize(ids::Type{<:IMAS.IDS}); kwargs...)
    txt = replace(split(split("$ids", "___")[end], "__")[end], r"\{.*\}" => "")
    return printstyled(io, txt; bold=true)
end

function AbstractTrees.printnode(io::IO, @nospecialize(ids::Type{<:IMAS.IDSvector}); kwargs...)
    txt = replace(split(split("$(eltype(ids))", "___")[end], "__")[end] * "[:]", r"\{.*\}" => "")
    return printstyled(io, txt; bold=true)
end

function AbstractTrees.printnode(io::IO, leaf::IMAS.IMASstructRepr; kwargs...)
    nfo = IMAS.info(leaf.ids_type, leaf.field)
    units = nfo.units
    if units == "-"
        units = ""
    else
        units = " [$(pretty_units(units))]"
    end
    return printstyled(io, "$(html_link_repr(leaf.field, leaf.location))$units")
end

function parameters_details_md(io::IO, pars::SimulationParameters.AbstractParameters)
    for leafRepr in AbstractTrees.Leaves(pars)
        leaf = leafRepr.value
        if typeof(leaf) <: SimulationParameters.ParametersVector
            error("$(SimulationParameters.spath(leaf)) has zero length, which prevents generation of documentation.")
        end
        if typeof(leaf) <: SimulationParameters.AbstractParameters
            continue
        end
        if ismissing(leaf.default)
            default = ""
            note = "note"
        else
            default = "* **Default:** `$(leaf.default)`"
            note = "tip"
        end
        if typeof(leaf) <: SimulationParameters.Switch
            options = "* **Options:** " * join(["`$(opt.first)`" for opt in leaf.options], ", ") * "\n    "
        else
            options = ""
        end
        txt = """

        ------------

        ```@raw html
        <div id='$(join(SimulationParameters.path(leaf),"."))'></div>
        ```
        !!! $note "$(join(SimulationParameters.path(leaf),"."))"
            $(leaf.description)
            * **Type:** `$(replace(string(typeof(leaf)),"SimulationParameters."=>""))`
            * **Units:** `$(isempty(leaf.units) ? "-" : leaf.units)`
            $(options)$(default)

        """
        write(io, txt)
    end
end

# ================== #
# generate main page #
# ================== #
include("$(@__DIR__)/src/index.jl")

# ================ #
# generate DD page #
# ================ #
include("$(@__DIR__)/src/dd_docs.jl")

# ==================== #
# generate Actors page #
# ==================== #
include("$(@__DIR__)/src/actors_docs.jl")

# =================== #
# generate inits page #
# =================== #
include("$(@__DIR__)/src/inits_docs.jl")

# ================= #
# generate ini page #
# ================= #
include("$(@__DIR__)/src/ini_docs.jl")

# ================= #
# generate act page #
# ================= #
include("$(@__DIR__)/src/act_docs.jl")

# =================== #
# generate cases page #
# =================== #
include("$(@__DIR__)/src/cases_docs.jl")

# ========================== #
# generate dependencies page #
# ========================== #
include("$(@__DIR__)/src/deps.jl")

# # ====================== #
# # generate examples page #
# # ====================== #
# include("$(@__DIR__)/src/examples.jl")

# ============== #
# build the docs #
# ============== #
makedocs(;
    root=@__DIR__,
    modules=[FUSE, IMAS, IMASdd],
    sitename="FUSE",
    build=joinpath(@__DIR__, "build"),
    format=Documenter.HTML(;
        repolink="https://github.com/ProjectTorreyPines/FUSE.jl",
        prettyurls=false,
        analytics="G-65D8V8C8VQ",
        sidebar_sitename=false,
        assets=["assets/favicon.ico"],
        size_threshold=nothing,
        size_threshold_warn=nothing
    ),
    repo=Remotes.GitHub("ProjectTorreyPines", "FUSE.jl"),
    warnonly=true,
    pages=[
        "Home" => "index.md",
        "Install" => "install.md",
        "Tutorial" => "tutorial.md",
        "Examples" => "examples.md",
        "Use Cases" => "cases.md",
        "Publications" => "pubs.md",
        "Actors" => "actors.md",
        "Initialization" => "inits.md",
        "`act` Parameters" => "act.md",
        "`ini` Parameters" => "ini.md",
        "`dd` Data Structure" => "dd.md",
        "Ecosystem" => "deps.md",
        "Development" => "develop.md",
        "License" => "license.md",
        "Notice" => "notice.md"
    ]
)

# convert "©(.*)©©(.*)©" patterns to hyperlinks
@info "Converting links"
for (file, parfile) in (("act", "act"), ("ini", "ini"), ("actors", "act"), ("dd", "dd"))
    local txt = open("$(@__DIR__)/build/$file.html", "r") do io
        return read(io, String)
    end
    txt = split(txt, "\n")
    for (k, line) in enumerate(txt)
        txt[k] = replace(txt[k], r"©(.*)©©©(.*)©" => s"<a href='©.html#\1'>\2</a>")
        txt[k] = replace(txt[k], r"©(.*)©©(.*)©" => s"<a href='©_details.html#\1'>\2</a>")
        txt[k] = replace(txt[k], "©" => parfile)
    end
    open("$(@__DIR__)/build/$file.html", "w") do io
        return write(io, join(txt, "\n"))
    end
end

# distinguish between input/output cells
@info "Styling examples"
for css in ("light", "dark")
    open("$(@__DIR__)/build/assets/themes/documenter-$css.css", "a") do io
        return write(
            io,
            """\n
.nohighlight {
background-color: transparent !important;
}"""
        )
    end
end
files_to_convert = readdir("$(@__DIR__)/build")[findall(x -> startswith(x, "example_") && endswith(x, ".html"), readdir("$(@__DIR__)/build"))]
for file in files_to_convert
    local txt = open("$(@__DIR__)/build/$file", "r") do io
        return read(io, String)
    end
    txt = split(txt, "\n")
    for (k, line) in enumerate(txt)
        txt[k] = replace(line, "<pre><code class=\"nohighlight" => "<pre class=\"nohighlight\"><code class=\"nohighlight")
    end
    open("$(@__DIR__)/build/$file", "w") do io
        return write(io, join(txt, "\n"))
    end
end

# Deploy docs
# This function deploys the documentation to the gh-pages branch of the repository.
# The main documentation that will be hosted on
# https://projecttorreypines.github.io/FUSE.jl/stable
# will be built from latest release tagged with a version number.
# The development documentation that will be hosted on
# https://projecttorreypines.github.io/FUSE.jl/dev
# will be built from the latest commit on the chosen devbranch argument below.
# For testing purposes, the devbranch argument can be set to WIP branch like "docs".
# While merging with master, the devbranch argument should be set to "master".
deploydocs(;
    repo="github.com/ProjectTorreyPines/FUSE.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "dev", "v#.#.#"]
)
