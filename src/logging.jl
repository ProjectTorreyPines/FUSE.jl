using Logging
import LoggingExtras

function uwanted_warnings_filter(log_args)
    return !any([
        contains(log_args.message, "both ImageMetadata and ImageAxes export"),
        startswith(log_args.message, "Keyword argument letter not supported with Plots.GRBackend")
    ])
end

const log_topics = Dict{Symbol,Logging.LogLevel}()
for topic in [:actors]
    log_topics[topic] = Logging.Error
end

"""
    logging(level::Logging.LogLevel, topic::Symbol, text::AbstractString)

Write to log with given level and topic
"""
function logging(level::Logging.LogLevel, topic::Symbol, text::AbstractString)
    if topic ∉ keys(log_topics)
        error("Valid FUSE logging topics are: $(Symbol[topic for topic in keys(log_topics)])")
    end
    Logging.@logmsg level text
end

"""
    logging(level::Union{Nothing,Logging.LogLevel}=nothing; topics...)

Set global ConsoleLogger based on level and set min_level for logging of different topics
"""
function logging(level::Union{Nothing,Logging.LogLevel}=nothing; topics...)
    if level !== nothing
        Logging.global_logger(LoggingExtras.ActiveFilteredLogger(uwanted_warnings_filter, Logging.ConsoleLogger(level)))
    end
    for (topic, value) in topics
        if topic ∈ keys(log_topics)
            if typeof(value) <: Logging.LogLevel
                log_topics[topic] = value
            else
                error("The value of a logging topic must be of type Logging.LogLevel")
            end
        else
            error("Valid FUSE logging topics are: $(Symbol[topic for topic in keys(log_topics)])")
        end
    end
    return nothing
end

function logging_actor_init(typeof_actor::DataType, args...; kw...)
    logging(Logging.Debug, :actors, "$typeof_actor @ init")
end

function step(actor::AbstractActor, args...; kw...)
    logging(Logging.Info, :actors, "$(typeof(actor)) @ step")
    return _step(actor, args...; kw...)
end

function finalize(actor::AbstractActor, args...; kw...)
    logging(Logging.Debug, :actors, "$(typeof(actor)) @finalize")
    return _finalize(actor, args...; kw...)
end
