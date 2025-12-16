#################################################
# Session state management for FUSE Web GUI     #
#################################################

using Logging
using Dates

###################
# Session Logger  #
###################

"""
    SessionLogger <: AbstractLogger

Custom logger that appends log messages to a Vector{String} buffer.
Formats messages with emoji indicators for different log levels.
"""
struct SessionLogger <: AbstractLogger
    buffer::Vector{String}
    min_level::LogLevel
end

SessionLogger(buffer::Vector{String}) = SessionLogger(buffer, Logging.Info)

Logging.min_enabled_level(logger::SessionLogger) = logger.min_level

Logging.shouldlog(logger::SessionLogger, level, _module, group, id) = level >= logger.min_level

function Logging.handle_message(logger::SessionLogger, level, message, _module, group, id,
                                file, line; kwargs...)
    level_str = if level == Logging.Debug
        "🔍 DEBUG"
    elseif level == Logging.Info
        "ℹ️ INFO"
    elseif level == Logging.Warn
        "⚠️ WARN"
    elseif level == Logging.Error
        "❌ ERROR"
    else
        "LOG"
    end

    msg = "$level_str: $message"
    push!(logger.buffer, msg)
end

Logging.catch_exceptions(logger::SessionLogger) = false


###################
# WebGuiSession   #
###################

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
    with_session_logging(f, session)

Execute function `f` while capturing all logging output, stdout, and stderr to the session's log buffer.
"""
function with_session_logging(f, session::WebGuiSession)
    logger = SessionLogger(session.log_buffer)

    # Create pipes for both stdout and stderr
    stdout_pipe = Pipe()
    stderr_pipe = Pipe()

    result = redirect_stdout(stdout_pipe) do
        redirect_stderr(stderr_pipe) do
            with_logger(logger) do
                f()
            end
        end
    end

    # Close write ends and read captured output
    close(Base.pipe_writer(stdout_pipe))
    close(Base.pipe_writer(stderr_pipe))

    stdout_output = String(read(stdout_pipe))
    stderr_output = String(read(stderr_pipe))

    close(Base.pipe_reader(stdout_pipe))
    close(Base.pipe_reader(stderr_pipe))

    # Add stdout output to log
    if !isempty(stdout_output)
        for line in split(stdout_output, '\n')
            if !isempty(strip(line))
                push!(session.log_buffer, strip(line))
            end
        end
    end

    # Add stderr output to log
    if !isempty(stderr_output)
        for line in split(stderr_output, '\n')
            if !isempty(strip(line))
                push!(session.log_buffer, strip(line))
            end
        end
    end

    return result
end

"""
    load_case!(session, case_sym, args, kwargs) -> Bool

Load a FUSE case into the session. Returns true on success, false on failure.
Logs progress and errors to session.log_buffer. Captures FUSE logger output.
"""
function load_case!(session::WebGuiSession, case_sym::Symbol, args::Tuple, kwargs::Dict{Symbol,Any})
    push!(session.log_buffer, "═══════════════════════════")
    push!(session.log_buffer, "Loading case: $case_sym")
    push!(session.log_buffer, "Initializing case parameters...")

    try
        # Create custom logger to capture @info/@warn/@error output
        logger = SessionLogger(session.log_buffer)

        # Create pipe to capture stderr
        pipe = Pipe()

        # Wrap FUSE calls with logger and stderr redirection
        redirect_stderr(pipe) do
            with_logger(logger) do
                session.ini, session.act = FUSE.case_parameters(Val{case_sym}(), args...; kwargs...)
                push!(session.log_buffer, "Creating IMAS data structure...")
                session.dd = IMAS.dd()

                push!(session.log_buffer, "Running FUSE.init!()...")
                FUSE.init!(session.dd, session.ini, session.act)
            end
        end

        # Close write end and read captured output
        close(Base.pipe_writer(pipe))
        stderr_output = String(read(pipe))
        close(Base.pipe_reader(pipe))

        # Add stderr output to log
        if !isempty(stderr_output)
            for line in split(stderr_output, '\n')
                if !isempty(strip(line))
                    push!(session.log_buffer, strip(line))
                end
            end
        end

        session.current_case = case_sym
        session.last_run_status = "success"
        push!(session.log_buffer, "✓ Successfully loaded $case_sym at $(Dates.now())")
        return true

    catch e
        bt = catch_backtrace()
        error_msg = sprint(showerror, e, bt)
        push!(session.log_buffer, "✗ Error loading $case_sym:")
        push!(session.log_buffer, error_msg)
        session.last_run_status = "error"
        return false
    end
end

"""
    is_loaded(session) -> Bool

Check if a case has been loaded into the session.
"""
is_loaded(session::WebGuiSession) = !isnothing(session.dd)
