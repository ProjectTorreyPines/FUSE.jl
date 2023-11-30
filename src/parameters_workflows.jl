abstract type ParametersFlow <: AbstractParameters end
abstract type AbstractWorkflow end

Base.@kwdef mutable struct ParametersFlowSettings{T} <: ParametersFlow where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :settings
    server::Switch{Symbol} = Switch{Symbol}([:localhost, :omega, :saga], "-", "Where to run"; default=:localhost)
end

Base.@kwdef mutable struct ParametersFlowTGLF_DB{T} <: ParametersFlow where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :FlowTGLF_DB
    settings::ParametersFlowSettings{T} = ParametersFlowSettings{T}()
    sat_rules::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "TGLF saturation rules to run")
    save_folder::Entry{String} = Entry{String}("-", "folder to save the database runs into"; default="")
    database_folder::Entry{String} = Entry{String}("-", "folder with input database"; default="")
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