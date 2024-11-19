using ProgressMeter
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
    dataframe::Union{DataFrame,Missing}
    iterator::Union{Vector{String},Missing}
    workflow::Union{Function,Missing}
end

function StudyDatabaseGenerator(sty::ParametersStudy, ini::ParametersAllInits, act::ParametersAllActors; kw...)
    sty = sty(kw...)
    study = StudyDatabaseGenerator(sty, ini, act, missing, missing, missing)
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

    # paraller run
    println("running $(sty.n_simulations) simulations with $(sty.n_workers) workers on $(sty.server)")
    results = @showprogress pmap(item -> run_case(study, item), iterator)

    # populate DataFrame
    for row in results
        if !isnothing(row)
            push!(study.dataframes_dict["outputs_summary"], row)
        end
    end

    analyze(study)

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
    extract_results(study)
    df = study.dataframe
    display(histogram(df.Te0; xlabel="Te0 [eV]", legend=false))
    display(histogram(df.Ti0; xlabel="Ti0 [eV]", legend=false))
    display(histogram(df.ne0; xlabel="ne0 [m⁻³]]", legend=false))
    return study
end

"""
    run_case(study::AbstractStudy, item::String)

Run a single case based by setting up a dd from ini and act and then executing the workflow_DatabaseGenerator workflow (feel free to change the workflow based on your needs)
"""
function run_case(study::AbstractStudy, item::String)
    act = study.act
    sty = study.sty
    @assert isa(study.workflow, Function) "Make sure to specicy a workflow to study.workflow that takes dd, ini , act as arguments"

    original_dir = pwd()
    savedir = abspath(joinpath(sty.save_folder, "$(item)__$(Dates.now())__$(getpid())"))
    mkdir(savedir)
    cd(savedir)

    # Redirect stdout and stderr to the file
    original_stdout = stdout  # Save the original stdout
    original_stderr = stderr  # Save the original stderr
    file_log = open("log.txt", "w")

    # deepcopy ini/act to avoid changes
    ini = rand(study.ini)
    dd = IMAS.dd()
    
    try
        redirect_stdout(file_log)
        redirect_stderr(file_log)

        study.workflow(dd, ini, act)

        # save simulation data to directory
        save(savedir, sty.save_dd ? dd : nothing, ini, act; timer=true, freeze=false, overwrite_files=true)

        return nothing
    catch error
        if isa(error, InterruptException)
            rethrow(error)
        end

        # save empty dd and error to directory
        save(savedir, nothing, ini, act; error, timer=true, freeze=false, overwrite_files=true)
    finally
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        cd(original_dir)
        close(file_log)
    end
end