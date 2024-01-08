abstract type ParametersFlow <: AbstractParameters end
abstract type AbstractWorkflow end

function flow_common_parameters(name::Symbol)
    if name == :server
        return Switch{String}(["localhost", "omega", "saga"], "-", "Where to run"; default="localhost")
    elseif name == :n_workers
        return Entry{Int}("-", "number of workers to run with")
    elseif name == :file_save_mode
        return Switch{Symbol}([:safe_write,:overwrite],"-", "The policy to implement when saving files, safe_write only writen when the folder is empty, overwrite overwrites"; default=:safe_write)

    end
end

Base.@kwdef mutable struct ParametersFlowTGLFdb{T} <: ParametersFlow where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :FlowTGLFdb
    server::Switch{String} = flow_common_parameters(:server)
    n_workers::Entry{Int} = flow_common_parameters(:n_workers)
    file_save_mode::Switch{Symbol} = flow_common_parameters(:file_save_mode)
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