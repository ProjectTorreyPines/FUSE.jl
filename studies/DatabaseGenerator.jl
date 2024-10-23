import JSON

#= ====================== =#
#  StudyDatabaseGenerator  #
#= ====================== =#

"""
    study_parameters(::Type{Val{:DatabaseGenerator}})::Tuple{FUSEparameters__ParametersStudyDatabaseGenerator,ParametersAllActors}

Generates a database of dds from ini and act based on ranges specified in ini (i.e. ini.equilibrium.R0 = 5.0 ↔ [4.0, 10.0])
There is a example notebook in FUSE_examples/study_database_generator.ipynb that goes through the steps of setting up, running and analyzing this study
"""
function study_parameters(::Type{Val{:DatabaseGenerator}})::Tuple{FUSEparameters__ParametersStudyDatabaseGenerator,ParametersAllActors}

    sty = FUSEparameters__ParametersStudyDatabaseGenerator{Real}()
    act = ParametersActors()

    # Change act for the default DatabaseGenerator run
    act.ActorCoreTransport.model = :FluxMatcher
    act.ActorFluxMatcher.evolve_pedestal = false
    act.ActorTGLF.warn_nn_train_bounds = false
    act.ActorFluxMatcher.evolve_rotation = :fixed


    # finalize 
    set_new_base!(sty)
    set_new_base!(act)

    return sty, act
end

Base.@kwdef mutable struct FUSEparameters__ParametersStudyDatabaseGenerator{T<:Real} <: ParametersStudy{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StudyDatabaseGenerator
    server::Switch{String} = study_common_parameters(; server="localhost")
    n_workers::Entry{Int} = study_common_parameters(; n_workers=missing)
    file_save_mode::Switch{Symbol} = study_common_parameters(; file_save_mode=:safe_write)
    release_workers_after_run::Entry{Bool} = study_common_parameters(; release_workers_after_run=true)
    save_dd::Entry{Bool} = study_common_parameters(; save_dd=true)
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the database runs into")
    n_simulations::Entry{Int} = Entry{Int}("-", "Number of sampled simulations")
end

mutable struct StudyDatabaseGenerator <: AbstractStudy
    sty::FUSEparameters__ParametersStudyDatabaseGenerator
    ini::ParametersAllInits
    act::ParametersAllActors
    dataframes_dict::Union{Dict{String,DataFrame},Missing}
    iterator::Union{Vector{String},Missing}
end

function StudyDatabaseGenerator(sty::ParametersStudy, ini::ParametersAllInits, act::ParametersAllActors; kw...)
    sty = sty(kw...)
    study = StudyDatabaseGenerator(sty, ini, act, missing, missing)
    return setup(study)
end

function _setup(study::StudyDatabaseGenerator)
    sty = study.sty

    check_and_create_file_save_mode(sty)

    parallel_environment(sty.server, sty.n_workers)

    return study
end

"""
    _run(study::StudyDatabaseGenerator)

Runs the DatabaseGenerator with sty settings in parallel on designated cluster
"""
function _run(study::StudyDatabaseGenerator)
    sty = study.sty

    @assert sty.n_workers == length(Distributed.workers()) "The number of workers =  $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"

    iterator = map(string, 1:sty.n_simulations)
    study.iterator = iterator
    study.dataframes_dict = Dict("outputs_summary" => StudyDatabaseGenerator_summary_dataframe() for name in study.iterator)

    # paraller run
    println("running $(sty.n_simulations) simulations with $(sty.n_workers) workers on $(sty.server)")
    results = pmap(item -> run_case(study, item), iterator)

    # populate DataFrame
    for row in results
        if !isnothing(row)
            push!(study.dataframes_dict["outputs_summary"], row)
        end
    end

    # Save JSON to a file
    json_data = JSON.json(study.dataframes_dict["outputs_summary"])
    open("$(sty.save_folder)/data_frame_outputs_summary.json", "w") do f
        return write(f, json_data)
    end

    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"
    end

    return study
end

"""
    _analyze(study::StudyDatabaseGenerator)

Example of analyze plots to display after the run feel free to change this method for your needs
"""
function _analyze(study::StudyDatabaseGenerator)
    display(histogram(study.dataframes_dict["outputs_summary"].Te0; xlabel="Te0 [eV]", legend=false))
    display(histogram(study.dataframes_dict["outputs_summary"].Ti0; xlabel="Ti0 [eV]", legend=false))
    display(histogram(study.dataframes_dict["outputs_summary"].ne0; xlabel="ne0 [m⁻³]]", legend=false))
    return study
end

"""
    run_case(study::AbstractStudy, item::String)

Run a single case based by setting up a dd from ini and act and then executing the workflow_DatabaseGenerator workflow (feel free to change the workflow based on your needs)
"""
function run_case(study::AbstractStudy, item::String)
    act = study.act
    sty = study.sty

    try
        # generate new ini
        ini = rand(study.ini)

        dd = IMAS.dd()
        workflow_DatabaseGenerator(dd, ini, act)

        if sty.save_dd
            IMAS.imas2json(dd, joinpath(sty.save_folder, "result_dd_$(item).json"))
            if parse(Bool, get(ENV, "FUSE_MEMTRACE", "false"))
                save(FUSE.memtrace, joinpath(sty.save_folder, "memtrace_$(item).txt"))
            end
        end

        return create_data_frame_row_DatabaseGenerator(dd)
    catch e
        if isa(e, InterruptException)
            rethrow(e)
        end
        open("$(sty.save_folder)/error_$(item).txt", "w") do file
            return showerror(file, e, catch_backtrace())
        end
    end
end

"""
    workflow_DatabaseGenerator(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Workflow of initalization and actors run on every dd to generate the database (feel free to change or write your own)
"""
function workflow_DatabaseGenerator(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    # initialize
    init(dd, ini, act)

    # Actors to run on the input dd
    actor_statplasma = ActorStationaryPlasma(dd, act)

    # whatever other actors you want to run can go here

    return actor_statplasma
end

"""
    create_data_frame_row_DatabaseGenerator(dd::IMAS.dd)

Example code of what to extract from the dd and store in the aggregate dataframe, FUSE.extract(dd) could be useful here as well.
All the dds will be stored in save_folder to run your own datagathering on
"""

function create_data_frame_row_DatabaseGenerator(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    return (ne0=cp1d.electrons.density_thermal[1], Te0=cp1d.electrons.temperature[1], Ti0=cp1d.t_i_average[1], zeff=cp1d.zeff[1])
end

function StudyDatabaseGenerator_summary_dataframe()
    return DataFrame(; ne0=Float64[], Te0=Float64[], Ti0=Float64[], zeff=Float64[])
end