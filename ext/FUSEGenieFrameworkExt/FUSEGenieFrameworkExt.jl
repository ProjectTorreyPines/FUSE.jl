"""
    FUSEGenieFrameworkExt

Julia package extension providing a web-based GUI for FUSE using GenieFramework.jl

This extension is automatically loaded when GenieFramework.jl is available in the environment.

# Main Entry Point
- `FUSE.launch_web_gui(; host="127.0.0.1", port=8050, open_browser=true)` - Launch the web GUI

# Dependencies
- GenieFramework.jl
  Meta package for Genie reactive apps. This packages exports
  Genie, Stipple, StippleUI, StipplePlotly, Stipple.Pages, Stipple.ModelStorage.Sessions, Stipple.ReactiveTools, Genie.Renderer.Html, Genie.Server
  and other packages from Genie Ecosystem as required in future.
"""
module FUSEGenieFrameworkExt

# Import main package
using FUSE
using IMAS
import FUSE: ParametersActors, ParametersInits
import SimulationParameters: Entry, Switch

# Import extension dependency
using GenieFramework

# Import standard library modules
using Dates

# Include submodules
include("session.jl")
include("parameters.jl")
include("plotting.jl")
include("ui.jl")

"""
    FUSE.launch_web_gui(; host="127.0.0.1", port=8050, open_browser=true)

Launch the FUSE web GUI using GenieFramework.jl.

# Keyword Arguments
- `host::String="127.0.0.1"` - Host address to bind the server to
- `port::Int=8050` - Port number for the web server
- `open_browser::Bool=true` - Automatically open the browser to the GUI

# Example
```julia
using Pkg
Pkg.add("GenieFramework")  # This also installs Stipple and Genie as dependencies

using FUSE
import GenieFramework
FUSE.launch_web_gui()  # Opens browser to http://127.0.0.1:8050
```

# Notes
- The server runs in blocking mode (async=false) by default
- Press Ctrl+C to stop the server
- Multiple users can connect to the same server instance
"""
function FUSE.launch_web_gui(;
    host::String="127.0.0.1",
    port::Int=8050,
    open_browser::Bool=true)

    # Register the main Stipple page
    @page("/", ui)

    # Display startup message
    url = "http://$host:$port"
    @info """
    ═══════════════════════════════════════════════════════════
    FUSE Web GUI starting...
    ═══════════════════════════════════════════════════════════
    URL: $url

    Press Ctrl+C to stop the server
    ═══════════════════════════════════════════════════════════
    """

    # Try to open browser if requested
    if open_browser
        try
            if Sys.islinux()
                run(`xdg-open $url`, wait=false)
            elseif Sys.isapple()
                run(`open $url`, wait=false)
            elseif Sys.iswindows()
                run(`cmd /c start $url`, wait=false)
            end
            @info "Opening browser..."
        catch e
            @warn "Could not open browser automatically" exception=e
            @info "Please open your browser to: $url"
        end
    else
        @info "Please open your browser to: $url"
    end

    # Start the Genie server (blocking)
    Genie.Server.up(host=host, port=port, async=false)
end

end # module FUSEGenieFrameworkExt
