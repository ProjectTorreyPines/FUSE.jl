using LoggingExtras

function uwanted_warnings_filter(log_args)
    return !any([
        contains(log_args.message, "both ImageMetadata and ImageAxes export"),
        startswith(log_args.message, "Keyword argument letter not supported with Plots.GRBackend")
        ])
end

logger = ActiveFilteredLogger(uwanted_warnings_filter, global_logger())
