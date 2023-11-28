abstract type ParametersFlow <: AbstractParameters end
abstract type AbstractWorkflow end

Base.@kwdef mutable struct ParametersFlowSettings <: ParametersFlow
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :settings
    server::Switch{Symbol} = Switch{Symbol}([:localhost, :omega, :saga], "-", "Where to run"; default=:localhost)
end

Base.@kwdef mutable struct ParametersFlowTGLF_DB{T} <: ParametersFlow where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :FlowTGLF_DB

    settings::ParametersFlowSettings=ParametersFlowSettings()
    ripple::Entry{T} = Entry{T}("-", "Fraction of toroidal field ripple evaluated at the outermost radius of the plasma chamber"; default=0.01)
    sat_rules::Entry{Union{AbstractVector{Symbol},Missing}} = Entry{Union{AbstractVector{Symbol},Missing}}("-", "sat rules to run"; default=missing)
    save_folder::Entry{String} = Entry{String}("-", "folder to save the database runs into"; default="")
    database_folder::Entry{String} = Entry{String}("-", "folder with input database"; default="")
end


#= ===== =#
#  flows  #
#= ===== =#

# NOTE only called once at precompile time, kernel needs to be restarted to include new file in cases
for filename in readdir(joinpath(@__DIR__, "..", "flows"))
    if endswith(filename, ".jl")
        @show filename
        include("../flows/" * filename)
    end
end


function flow_parameters(flow::Symbol; kw...)
    if length(methods(flow_parameters, (Type{Val{flow}},))) == 0
        throw(InexistentParameterException([flow]))
    end
    return flow_parameters(Val{flow})#; kw...)
end
