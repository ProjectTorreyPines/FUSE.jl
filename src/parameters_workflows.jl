abstract type ParametersFlow <: AbstractParameters end
abstract type AbstractWorkflow end

"""
    flow_common_parameters(; kw...)

Returns commonly used parameters as a switch or entry, example: flow_common_parameters(server="localhost")
"""
function flow_common_parameters(; kw...)
    @assert length(kw) == 1 "flow_common_parameters only takes one argument"
    name = first(keys(kw))
    default = first(values(kw))
    if name == :server
        return Switch{String}(["localhost", "omega", "saga"], "-", "Where to run"; default)
    elseif name == :n_workers
        return Entry{Int}("-", "number of workers to run with"; default)
    elseif name == :file_save_mode
        return Switch{Symbol}([:safe_write, :overwrite], "-", "Saving file policy, safe_write only writes when the folder is empty"; default)
    elseif name == :release_workers_after_run
        return Entry{Bool}("-", "releases the workers after running the workflow"; default)
    elseif name == :keep_output_dd
        return Entry{Bool}("-", "Store the output dds of the workflow run"; default)
    else
        error("There is no flow_common_parameter for name = $name")
    end
end

Base.@kwdef mutable struct ParametersFlowTGLFdb{T} <: ParametersFlow where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :FlowTGLFdb
    server::Switch{String} = flow_common_parameters(; server="localhost")
    n_workers::Entry{Int} = flow_common_parameters(; n_workers=missing)
    file_save_mode::Switch{Symbol} = flow_common_parameters(; file_save_mode=:safe_write)
    release_workers_after_run::Entry{Bool} = flow_common_parameters(; release_workers_after_run=true)
    keep_output_dd::Entry{Bool} = flow_common_parameters(; keep_output_dd=true)
    sat_rules::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "TGLF saturation rules to run")
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the database runs into")

    database_folder::Entry{String} = Entry{String}("-", "Folder with input database")
end


#= ===== =#
#  flows  #
#= ===== =#

# NOTE only called once at precompile time, kernel needs to be restarted to include new file in `flows` directory
for filename in readdir(joinpath(@__DIR__, "..", "flows"))
    if endswith(filename, ".jl")
        include("../flows/" * filename)
    end
end

"""
    function name(wf::AbstractWorkflow)

Returns the name of the workflow as a string
"""
function name(wf::AbstractWorkflow)
    return string(split(string(typeof(wf)), ".")[end])
end

function flow_parameters(flow::Symbol; kw...)
    if length(methods(flow_parameters, (Type{Val{flow}},))) == 0
        throw(InexistentParameterException([flow]))
    end
    return flow_parameters(Val{flow}; kw...)
end

function setup(wf::AbstractWorkflow)
    return _setup(wf)
end

function analyze(wf::AbstractWorkflow)
    # here you can add timing info and more
    return _analyze(wf)
end

function run(wf::T, args...; kw...) where {T<:AbstractWorkflow}
    timer_name = name(wf)
    TimerOutputs.reset_timer!(timer_name)
    TimerOutputs.@timeit timer timer_name begin
        TimerOutputs.reset_timer!(timer_name)
        memory_time_tag("workflow $(name(wf)) - @start")
        wf = _run(wf, args...; kw...)
        memory_time_tag("workflow $(name(wf)) - @finish")
    end
    return wf
end


"""
    check_and_create_file_save_mode(flw)

Checks the selected file_save_mode and creates the folder accordingly
"""
function check_and_create_file_save_mode(flw)
    @assert !ismissing(getproperty(flw, :save_folder, missing)) "Make sure flw.save_folder = $(flw.save_folder) is set"
    if flw.file_save_mode == :safe_write
        if isdir(flw.save_folder)
            @assert isempty(readdir((flw.save_folder))) "$(flw.save_folder) isn't empty, change flw.file_save_mode or point to a new save folder"
        else
            @assert !isfile(flw.save_folder) "$(flw.save_folder) can't be a file"
            mkdir(flw.save_folder)
        end
    elseif flw.file_save_mode == :overwrite && !isdir(flw.save_folder)
        mkdir(flw.save_folder)
    end
end