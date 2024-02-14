abstract type ParametersStudy{T} <: AbstractParameters{T} end
abstract type AbstractStudy end

"""
    study_common_parameters(; kw...)

Returns commonly used parameters as a switch or entry, example: study_common_parameters(server="localhost")
"""
function study_common_parameters(; kw...)
    @assert length(kw) == 1 "study_common_parameters only takes one argument"
    name = first(keys(kw))
    default = first(values(kw))
    if name == :server
        return Switch{String}(["localhost", "omega", "saga"], "-", "Where to run"; default)
    elseif name == :n_workers
        return Entry{Int}("-", "number of workers to run with"; default)
    elseif name == :file_save_mode
        return Switch{Symbol}([:safe_write, :overwrite], "-", "Saving file policy, safe_write only writes when the folder is empty"; default)
    elseif name == :release_workers_after_run
        return Entry{Bool}("-", "releases the workers after running the study"; default)
    elseif name == :keep_output_dd
        return Entry{Bool}("-", "Store the output dds of the study run"; default)
    else
        error("There is no study_common_parameter for name = $name")
    end
end

#= ============ =#
#  studies  #
#= ============ =#

# NOTE only called once at precompile time, kernel needs to be restarted to include new file in `studies` directory
for filename in readdir(joinpath(@__DIR__, "..", "studies"))
    if endswith(filename, ".jl")
        include("../studies/" * filename)
    end
end

"""
    function name(study::AbstractStudy)

Returns the name of the study as a string
"""
function name(study::AbstractStudy)
    return string(split(string(typeof(study)), ".")[end])
end

function study_parameters(study::Symbol; kw...)
    if length(methods(study_parameters, (Type{Val{study}},))) == 0
        throw(InexistentParameterException([study]))
    end
    return study_parameters(Val{study}; kw...)
end

function setup(study::AbstractStudy)
    return _setup(study)
end

function analyze(study::AbstractStudy)
    # here you can add timing info and more
    return _analyze(study)
end

function run(study::T, args...; kw...) where {T<:AbstractStudy}
    timer_name = name(study)
    TimerOutputs.reset_timer!(timer_name)
    TimerOutputs.@timeit timer timer_name begin
        TimerOutputs.reset_timer!(timer_name)
        memory_time_tag("study $(name(study)) - @start")
        wf = _run(study, args...; kw...)
        memory_time_tag("study $(name(study)) - @finish")
    end
    return wf
end


"""
    check_and_create_file_save_mode(sty)

Checks the selected file_save_mode and creates the folder accordingly
"""
function check_and_create_file_save_mode(sty)
    @assert !ismissing(getproperty(sty, :save_folder, missing)) "Make sure sty.save_folder = $(sty.save_folder) is set"
    if sty.file_save_mode == :safe_write
        if isdir(sty.save_folder)
            @assert isempty(readdir((sty.save_folder))) "$(sty.save_folder) isn't empty, change sty.file_save_mode or point to a new save folder"
        else
            @assert !isfile(sty.save_folder) "$(sty.save_folder) can't be a file"
            mkdir(sty.save_folder)
        end
    elseif sty.file_save_mode == :overwrite && !isdir(sty.save_folder)
        mkdir(sty.save_folder)
    end
end