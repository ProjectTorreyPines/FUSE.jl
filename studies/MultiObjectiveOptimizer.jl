#= ============================ =#
#  StudyMultiObjectiveOptimizer  #
#= ============================ =#

"""
    study_parameters(::Type{Val{:MultiObjectiveOptimizer}})::Tuple{FUSEparameters__ParametersStudyMultiObjectiveOptimizer,ParametersAllActors}

Generates a database of dds from ini and act based on ranges specified in ini
"""
function study_parameters(::Type{Val{:MultiObjectiveOptimizer}})::Tuple{FUSEparameters__ParametersStudyMultiObjectiveOptimizer,ParametersAllActors}

    sty = FUSEparameters__ParametersStudyMultiObjectiveOptimizer{Real}()
    act = ParametersActors()

    # finalize 
    set_new_base!(sty)
    set_new_base!(act)

    return sty, act
end

Base.@kwdef mutable struct FUSEparameters__ParametersStudyMultiObjectiveOptimizer{T<:Real} <: ParametersStudy{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StudyMultiObjectiveOptimizer
    # HPC related parameters
    server::Switch{String} = study_common_parameters(; server="localhost")
    n_workers::Entry{Int} = study_common_parameters(; n_workers=missing)
    file_save_mode::Switch{Symbol} = study_common_parameters(; file_save_mode=:safe_write)
    release_workers_after_run::Entry{Bool} = study_common_parameters(; release_workers_after_run=true)
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the database runs into")
    # Optimization related parameters
    population_size::Entry{Int} = Entry{Int}("-", "Number of individuals in a generation")
    number_of_generations::Entry{Int} = Entry{Int}("-", "Number generations")
#    worfklow:: = Entry{}("", default=ActorWholeFacility)
end

mutable struct StudyMultiObjectiveOptimizer <: AbstractStudy
    sty::FUSEparameters__ParametersStudyMultiObjectiveOptimizer
    ini::ParametersAllInits
    act::ParametersAllActors
    constraint_functions::Vector{IMAS.ConstraintFunction}
    objective_functions::Vector{IMAS.ObjectiveFunction}
    dataframe::Union{DataFrame,Missing}
end

function StudyMultiObjectiveOptimizer(sty::ParametersStudy, ini::ParametersAllInits, act::ParametersAllActors, constraint_functions::Vector{IMAS.ConstraintFunction}, objective_functions::Vector{IMAS.ObjectiveFunction}; kw...)
    sty = sty(kw...)
    study = StudyMultiObjectiveOptimizer(sty, ini, act, constraint_functions, objective_functions, missing)
    return setup(study)
end

function _setup(study::StudyMultiObjectiveOptimizer)
    sty = study.sty

    check_and_create_file_save_mode(sty)

    parallel_environment(sty.server, sty.n_workers)

    return study
end

function _run(study::StudyMultiObjectiveOptimizer)
    sty = study.sty

    @assert sty.n_workers == length(Distributed.workers()) "The number of workers =  $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"
    @assert iseven(sty.number_of_generations) "Population size must be even"

    optimization_parameters = Dict(
        :N => sty.population_size,
        :iterations => sty.number_of_generations,
        :continue_state => nothing,
        :save_folder => sty.save_folder)
    @show optimization_parameters
    state = workflow_multiobjective_optimization(
        study.ini, study.act, ActorWholeFacility, study.objective_functions,
        study.constraint_functions; optimization_parameters...)

    save_optimization(
        joinpath(sty.save_folder, "optimization_state.bson"),
        state,
        study.ini,
        study.act,
        study.objective_functions,
        study.constraint_functions)

    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"
    end

    return study
end

function _analyze(study::StudyMultiObjectiveOptimizer)
    prinln("analyzing study")
    return study
end