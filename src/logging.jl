using Logging
import LoggingExtras

function uwanted_warnings_filter(log_args)
    return !any([
        contains(log_args.message, "both ImageMetadata and ImageAxes export"),
        startswith(log_args.message, "Keyword argument letter not supported with Plots.GRBackend")
    ])
end

logger = LoggingExtras.ActiveFilteredLogger(uwanted_warnings_filter, Logging.global_logger())

mutable struct LogOptions
    actors::Bool
end

const logging_options = LogOptions(false)

function logging(typeof_actor::DataType, args...; kw...)
    #logging_text(typeof_actor, "INIT")
    return nothing
end

function step(actor::AbstractActor, args...; kw...)
    #logging_text(actor, "STEP")
    if logging_options.actors
        logging_text(actor)
    end
    return _step(actor, args...; kw...)
end

function finalize(actor::AbstractActor, args...; kw...)
    #logging_text(actor, "FINALIZE")
    return _finalize(actor, args...; kw...)
end

function logging_text(actor::AbstractActor, text="")
    return logging_text(typeof(actor), text)
end

function logging_text(typeof_actor::DataType, text="")
    actor_name = replace("$typeof_actor", "FUSE.Actor" => "")
    if !isempty(text)
        text = ": $text"
    end
    @info "$actor_name$text"
    return nothing
end