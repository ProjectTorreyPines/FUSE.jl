module WeaveExt

using FUSE
using FUSE: IMAS, ParametersAllInits, ParametersAllActors, WeaveBackend
import Weave
using Logging

"""
    weave_digest(::WeaveBackend, dd::IMAS.dd, title::AbstractString, description::AbstractString=""; 
                ini::Union{Nothing,ParametersAllInits}=nothing,
                act::Union{Nothing,ParametersAllActors}=nothing)

Generate PDF digest using Weave. This method is only available when Weave.jl is loaded.

# Arguments
- `::WeaveBackend`: Dispatch type indicating Weave functionality is available
- `dd::IMAS.dd`: IMAS data structure containing simulation results
- `title::AbstractString`: Title for the PDF report
- `description::AbstractString=""`: Optional description for the report
- `ini::Union{Nothing,ParametersAllInits}=nothing`: Optional initialization parameters
- `act::Union{Nothing,ParametersAllActors}=nothing`: Optional actor parameters

# Returns
- `AbstractString`: Path to the generated PDF file, or `nothing` if generation failed

# Notes
This function processes the Weave template files (`digest.jmd` and `digest.tpl`) located
in the same directory as this extension to generate a comprehensive PDF report.
"""
function FUSE.weave_digest(::WeaveBackend, dd::IMAS.dd,
    title::AbstractString,
    description::AbstractString="";
    ini::Union{Nothing,ParametersAllInits}=nothing,
    act::Union{Nothing,ParametersAllActors}=nothing
)
    title = replace(title, r".pdf$" => "", "_" => " ")
    outfilename = joinpath(pwd(), "$(replace(title," "=>"_")).pdf")

    tmpdir = mktempdir()
    logger = SimpleLogger(stderr, Logging.Warn)
    try
        filename = redirect_stdout(Base.DevNull()) do
            filename = with_logger(logger) do
                return Weave.weave(joinpath(@__DIR__, "digest.jmd");
                    latex_cmd=["xelatex"],
                    mod=FUSE,
                    doctype="md2pdf",
                    template=joinpath(@__DIR__, "digest.tpl"),
                    out_path=tmpdir,
                    args=Dict(
                        :dd => dd,
                        :ini => ini,
                        :act => act,
                        :title => title,
                        :description => description))
            end
        end
        cp(filename, outfilename; force=true)
        return outfilename
    catch e
        if isa(e, InterruptException)
            rethrow(e)
        end
        println("Generation of $(basename(outfilename)) failed. See directory: $tmpdir\n$e")
        return nothing
    else
        rm(tmpdir; recursive=true, force=true)
    end
end

end # module WeaveExt