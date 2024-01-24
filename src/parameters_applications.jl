abstract type ParametersApplication <: AbstractParameters end
abstract type AbstractApplication end

"""
    application_common_parameters(; kw...)

Returns commonly used parameters as a switch or entry, example: application_common_parameters(server="localhost")
"""
function application_common_parameters(; kw...)
    @assert length(kw) == 1 "application_common_parameters only takes one argument"
    name = first(keys(kw))
    default = first(values(kw))
    if name == :server
        return Switch{String}(["localhost", "omega", "saga"], "-", "Where to run"; default)
    elseif name == :n_workers
        return Entry{Int}("-", "number of workers to run with"; default)
    elseif name == :file_save_mode
        return Switch{Symbol}([:safe_write, :overwrite], "-", "Saving file policy, safe_write only writes when the folder is empty"; default)
    elseif name == :release_workers_after_run
        return Entry{Bool}("-", "releases the workers after running the application"; default)
    elseif name == :keep_output_dd
        return Entry{Bool}("-", "Store the output dds of the application run"; default)
    else
        error("There is no application_common_parameter for name = $name")
    end
end

#= ============ =#
#  applications  #
#= ============ =#

# NOTE only called once at precompile time, kernel needs to be restarted to include new file in `applications` directory
for filename in readdir(joinpath(@__DIR__, "..", "applications"))
    if endswith(filename, ".jl")
        include("../applications/" * filename)
    end
end

"""
    function name(application::AbstractApplication)

Returns the name of the application as a string
"""
function name(application::AbstractApplication)
    return string(split(string(typeof(application)), ".")[end])
end

function application_parameters(application::Symbol; kw...)
    if length(methods(application_parameters, (Type{Val{application}},))) == 0
        throw(InexistentParameterException([application]))
    end
    return application_parameters(Val{application}; kw...)
end

function setup(application::AbstractApplication)
    return _setup(application)
end

function analyze(application::AbstractApplication)
    # here you can add timing info and more
    return _analyze(application)
end

function run(application::T, args...; kw...) where {T<:AbstractApplication}
    timer_name = name(application)
    TimerOutputs.reset_timer!(timer_name)
    TimerOutputs.@timeit timer timer_name begin
        TimerOutputs.reset_timer!(timer_name)
        memory_time_tag("application $(name(application)) - @start")
        wf = _run(application, args...; kw...)
        memory_time_tag("application $(name(application)) - @finish")
    end
    return wf
end


"""
    check_and_create_file_save_mode(app)

Checks the selected file_save_mode and creates the folder accordingly
"""
function check_and_create_file_save_mode(app)
    @assert !ismissing(getproperty(app, :save_folder, missing)) "Make sure app.save_folder = $(app.save_folder) is set"
    if app.file_save_mode == :safe_write
        if isdir(app.save_folder)
            @assert isempty(readdir((app.save_folder))) "$(app.save_folder) isn't empty, change app.file_save_mode or point to a new save folder"
        else
            @assert !isfile(app.save_folder) "$(app.save_folder) can't be a file"
            mkdir(app.save_folder)
        end
    elseif app.file_save_mode == :overwrite && !isdir(app.save_folder)
        mkdir(app.save_folder)
    end
end