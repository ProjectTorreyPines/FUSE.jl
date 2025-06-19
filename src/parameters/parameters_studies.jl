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
