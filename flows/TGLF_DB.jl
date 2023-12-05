"""
    flow_parameters(::Type{Val{:TGLF_DB}})::Tuple{ParametersFlowTGLF_DB, ParametersAllActors}
"""
function flow_parameters(::Type{Val{:TGLF_DB}})::Tuple{ParametersFlowTGLF_DB,ParametersAllActors}

    flw = ParametersFlowTGLF_DB{Real}()
    act = ParametersActors()

    #### FLW ####
    flw.settings.server = :localhost

    # finalize 
    set_new_base!(flw)
    set_new_base!(act)

    return flw, act
end

mutable struct FlowTGLF_DB <: AbstractWorkflow
    flw::ParametersFlowTGLF_DB
    act::ParametersAllActors
    dataframe::DataFrame
end

function FlowTGLF_DB(flw, act)
    return FlowTGLF_DB(flw, act, DataFrames.DataFrame())
end

function _setup(wf::FlowTGLF_DB)
    println("setting up")
    return wf
end

function _run(wf::FlowTGLF_DB)
    println("running FlowTGLF_DB")
    return wf
end

function _analyze(wf::FlowTGLF_DB)
    println("analyzing FlowTGLF_DB")
    return wf
end

