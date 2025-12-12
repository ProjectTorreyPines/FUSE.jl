#################################################
# Session state management for FUSE Web GUI     #
#################################################

"""
    WebGuiSession

Mutable struct holding the session state for a web GUI user.
Stored as a @private variable so it's not synced to the browser.
"""
mutable struct WebGuiSession
    current_case::Union{Nothing,Symbol}
    dd::Union{Nothing,IMAS.dd}
    ini::Union{Nothing,FUSE.ParametersInits}
    act::Union{Nothing,FUSE.ParametersActors}
    log_buffer::Vector{String}
    last_run_status::String
end

WebGuiSession() = WebGuiSession(nothing, nothing, nothing, nothing, String[], "")

"""
    append_log!(session, msg)

Add a message to the session's log buffer.
"""
function append_log!(session::WebGuiSession, msg::AbstractString)
    push!(session.log_buffer, msg)
    return nothing
end

"""
    clear_log!(session)

Clear the session's log buffer.
"""
function clear_log!(session::WebGuiSession)
    empty!(session.log_buffer)
    return nothing
end

"""
    load_case!(session, case_sym, args, kwargs) -> Bool

Load a FUSE case into the session. Returns true on success, false on failure.
Logs progress and errors to session.log_buffer.
"""
function load_case!(session::WebGuiSession, case_sym::Symbol, args::Tuple, kwargs::Dict{Symbol,Any})
    append_log!(session, "Loading case: $case_sym")

    try
        # Call case_parameters with the case Val type
        session.ini, session.act = FUSE.case_parameters(Val{case_sym}(), args...; kwargs...)
        append_log!(session, "✓ case_parameters completed")

        # Create and initialize dd
        session.dd = IMAS.dd()
        FUSE.init!(session.dd, session.ini, session.act)
        append_log!(session, "✓ init! completed")

        session.current_case = case_sym
        session.last_run_status = "success"
        return true

    catch e
        error_msg = sprint(showerror, e)
        append_log!(session, "✗ Error: $error_msg")
        session.last_run_status = "error"
        return false
    end
end

"""
    is_loaded(session) -> Bool

Check if a case has been loaded into the session.
"""
is_loaded(session::WebGuiSession) = !isnothing(session.dd)
