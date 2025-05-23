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
        return Switch{String}(["localhost", "omega", "saga", "feynman", "engaging"], "-", "Where to run"; default)
    elseif name == :n_workers
        return Entry{Int}("-", "Number of workers to run with"; default)
    elseif name == :file_save_mode
        return Switch{Symbol}([:safe_write, :overwrite, :append], "-", "Saving file policy, `safe_write` only writes when the folder is empty"; default)
    elseif name == :release_workers_after_run
        return Entry{Bool}("-", "Releases the workers after running the study"; default)
    elseif name == :save_dd
        return Entry{Bool}("-", "Save dd of the study to save folder"; default)
    elseif name == :database_policy
        return Switch{Symbol}([:separate_folders, :single_hdf5], "-", "Data storage policy: 'separate_folders' stores each case in a separate folder, while 'single_hdf5' merges all cases into a single HDF5 file"; default)
    else
        error("There is no study_common_parameter named `$name`")
    end
end

#= ======= =#
#  studies  #
#= ======= =#

# NOTE only called once at precompile time, kernel needs to be restarted to include new file in `studies` directory
for filename in readdir(joinpath(@__DIR__, "..", "..", "studies"))
    if endswith(filename, ".jl")
        include(joinpath(@__DIR__, "..", "..", "studies", filename))
    end
end

"""
    name(study::AbstractStudy)

Returns the name of the study as a string
"""
function name(study::AbstractStudy)
    return string(split(string(typeof(study)), ".")[end])
end

function study_parameters(study::Symbol; kw...)
    if length(methods(study_parameters, (Type{Val{study}},))) == 0
        error("study `$study` does not exist.\nPossible options are:\n\n$(join(["$method" for method in methods(study_parameters)],"\n"))")
    end
    return study_parameters(Val{study}; kw...)
end

function setup(study::AbstractStudy)
    return _setup(study)
end

function analyze(study::AbstractStudy; kw...)
    # here you can add timing info and more
    return _analyze(study; kw...)
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
    @assert !isempty(sty.save_folder) "Make sure sty.save_folder = $(sty.save_folder) is set"
    if sty.file_save_mode == :safe_write
        if isdir(sty.save_folder)
            @assert isempty(readdir((sty.save_folder))) "$(sty.save_folder) isn't empty, change sty.file_save_mode or point to a new save folder"
        else
            @assert !isfile(sty.save_folder) "$(sty.save_folder) can't be a file"
            mkdir(sty.save_folder)
        end
    elseif sty.file_save_mode == :overwrite
        rm(sty.save_folder; force=true, recursive=true)
        mkdir(sty.save_folder)
    elseif sty.file_save_mode == :append
        # this is fine
    end
end
