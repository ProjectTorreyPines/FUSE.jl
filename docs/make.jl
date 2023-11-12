using Pkg
Pkg.activate(".")
using Documenter
import FUSE
import IMAS
import IMASDD
import SimulationParameters
import AbstractTrees
import ProgressMeter
import Dates
using InteractiveUtils: subtypes

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
    info = IMAS.info(leaf.location)
    units = get(info, "units", "-")
    if units == "-"
        units = ""
    else
        units = " [$(pretty_units(units))]"
    end
    return printstyled(io, "$(html_link_repr(leaf.key, leaf.location))$units")
end

function parameters_details_md(io::IO, pars::SimulationParameters.AbstractParameters)
    for leafRepr in AbstractTrees.Leaves(pars)
        leaf = leafRepr.value
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
        if typeof(leaf) <: FUSE.Switch
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

# ====================== #
# generate examples page #
# ====================== #
include("$(@__DIR__)/src/examples.jl")

# ============== #
# build the docs #
# ============== #
makedocs(;
    root=@__DIR__,
    modules=[FUSE, IMAS, IMASDD],
    sitename="FUSE",
    format=Documenter.HTML(; prettyurls=false, sidebar_sitename=false, assets=["assets/favicon.ico"]),
    remotes=nothing,
    warnonly=true,
    pages=[
        "Concepts" => "index.md",
        "Data Structure" => "dd.md",
        "Actors" => "actors.md",
        "Parameters" => ["ini Parameters" => "ini.md", "act Parameters" => "act.md", "Use Cases" => "cases.md", "Initialization" => "inits.md"],
        "Examples" => "examples.md",
        "Development" => "develop.md",
        "Install" => ["Install FUSE" => "install.md", "on SAGA" => "install_saga.md", "on OMEGA" => "install_omega.md"],
        "Others" => ["GASC" => "gasc.md", "Utilities" => "utils.md"]
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
