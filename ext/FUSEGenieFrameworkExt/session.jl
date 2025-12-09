# Session state management for FUSE web GUI

using Logging

# Custom logger that appends to session log buffer
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

function Logging.min_enabled_level(logger::SessionLogger)
    return logger.min_level
end

function Logging.shouldlog(logger::SessionLogger, level, _module, group, id)
    return level >= logger.min_level
end

function Logging.handle_message(logger::SessionLogger, level, message, _module, group, id,
                                 file, line; kwargs...)
    # Format log message with emoji indicators
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

"""
    WebGuiSession

Mutable struct holding the state of a FUSE web GUI session.

# Fields
- `current_case::Union{Nothing, Symbol}` - Currently loaded case key
- `case_args::Dict{String, Any}` - Case positional arguments
- `case_kwargs::Dict{Symbol, Any}` - Case keyword arguments
- `dd::Union{Nothing, IMAS.dd}` - IMAS data structure
- `ini::Union{Nothing, ParametersInits}` - Initialization parameters
- `act::Union{Nothing, ParametersActors}` - Actor parameters
- `selected_actor::Union{Nothing, Symbol}` - Currently selected actor
- `last_run_status::String` - Status of last actor run ("success", "error", or "")
- `last_run_time::Union{Nothing, DateTime}` - Timestamp of last actor run
- `log_buffer::Vector{String}` - Log messages
- `case_cache::Dict` - Cache for loaded cases
"""
mutable struct WebGuiSession
    # Current case
    current_case::Union{Nothing,Symbol}
    case_args::Dict{String,Any}
    case_kwargs::Dict{Symbol,Any}

    # FUSE data structures
    dd::Union{Nothing,IMAS.dd}
    ini::Union{Nothing,ParametersInits}
    act::Union{Nothing,ParametersActors}

    # Actor execution
    selected_actor::Union{Nothing,Symbol}
    last_run_status::String
    last_run_time::Union{Nothing,DateTime}

    # Logging
    log_buffer::Vector{String}

    # Cache for loaded cases
    case_cache::Dict{Tuple{Symbol,Any,Any},Tuple{Any,Any,Any}}
end

"""
    WebGuiSession()

Create a new WebGuiSession with empty/default values.
"""
function WebGuiSession()
    WebGuiSession(
        nothing, Dict(), Dict(),
        nothing, nothing, nothing,
        nothing, "", nothing,
        String[],
        Dict()
    )
end

"""
    load_case!(session::WebGuiSession, case_sym::Symbol, args, kwargs)

Load a FUSE case into the session. Uses caching to avoid recompilation.

# Arguments
- `session::WebGuiSession` - The session to modify
- `case_sym::Symbol` - The case key (e.g., :KDEMO, :ITER)
- `args` - Positional arguments for case_parameters
- `kwargs` - Keyword arguments for case_parameters

# Returns
- `Bool` - true if successful, false otherwise
"""
function load_case!(session::WebGuiSession, case_sym::Symbol, args, kwargs)
    # Add separator and starting message
    push!(session.log_buffer, "═══════════════════════════")
    push!(session.log_buffer, "Loading case: $case_sym")

    cache_key = (case_sym, args, kwargs)

    # Check cache first
    if haskey(session.case_cache, cache_key)
        session.dd, session.ini, session.act = deepcopy(session.case_cache[cache_key])
        session.current_case = case_sym
        session.case_args = args
        session.case_kwargs = kwargs
        push!(session.log_buffer, "✓ Reloaded from cache at $(Dates.now())")
        return true
    end

    # Load fresh case
    push!(session.log_buffer, "Initializing case parameters...")

    try
        # Create custom logger to capture @info/@warn/@error output
        logger = SessionLogger(session.log_buffer)

        # Create pipe to capture stderr (FUSE's custom logging)
        pipe = Pipe()

        # Wrap FUSE calls with logger and stderr redirection
        redirect_stderr(pipe) do
            with_logger(logger) do
                session.ini, session.act = FUSE.case_parameters(case_sym; kwargs...)
                push!(session.log_buffer, "Creating IMAS data structure...")
                session.dd = IMAS.dd()

                push!(session.log_buffer, "Running FUSE.init()...")
                FUSE.init(session.dd, session.ini, session.act)
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

        # Cache the loaded case
        session.case_cache[cache_key] = (session.dd, session.ini, session.act)
        session.current_case = case_sym
        session.case_args = args
        session.case_kwargs = kwargs

        push!(session.log_buffer, "✓ Successfully loaded $case_sym at $(Dates.now())")
        return true
    catch e
        bt = catch_backtrace()
        error_msg = sprint(showerror, e, bt)
        push!(session.log_buffer, "✗ Error loading $case_sym:")
        push!(session.log_buffer, error_msg)
        return false
    end
end

"""
    run_actor!(session::WebGuiSession, actor_sym::Symbol)

Run a FUSE actor on the current session's data structure.

# Arguments
- `session::WebGuiSession` - The session containing dd, act
- `actor_sym::Symbol` - The actor to run (e.g., :ActorFluxMatcher)

# Returns
- `Bool` - true if successful, false otherwise
"""
function run_actor!(session::WebGuiSession, actor_sym::Symbol)
    # Check prerequisites
    if isnothing(session.dd) || isnothing(session.act)
        push!(session.log_buffer, "✗ No case loaded. Please load a case first.")
        session.last_run_status = "error"
        return false
    end

    # Add separator and starting message
    push!(session.log_buffer, "═══════════════════════════")
    push!(session.log_buffer, "Running actor: $actor_sym")

    try
        # Create custom logger to capture @info/@warn/@error output
        logger = SessionLogger(session.log_buffer)

        # Create pipe to capture stderr (FUSE's custom logging)
        pipe = Pipe()

        # Wrap actor calls with logger and stderr redirection
        redirect_stderr(pipe) do
            with_logger(logger) do
                if actor_sym == :ActorFluxMatcher
                    FUSE.ActorFluxMatcher(session.dd, session.act)
                elseif actor_sym == :ActorNoOperation
                    FUSE.ActorNoOperation(session.dd, session.act)
                else
                    error("Unknown actor: $actor_sym. Supported actors: :ActorNoOperation, :ActorFluxMatcher")
                end
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

        # Update session state
        session.selected_actor = actor_sym
        session.last_run_time = Dates.now()
        session.last_run_status = "success"
        push!(session.log_buffer, "✓ Successfully ran $actor_sym at $(session.last_run_time)")
        return true
    catch e
        bt = catch_backtrace()
        error_msg = sprint(showerror, e, bt)
        session.last_run_status = "error"
        push!(session.log_buffer, "✗ Error running $actor_sym:")
        push!(session.log_buffer, error_msg)
        return false
    end
end
