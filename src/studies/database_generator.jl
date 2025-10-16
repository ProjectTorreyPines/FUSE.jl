import HDF5

#= ====================== =#
#  StudyDatabaseGenerator  #
#= ====================== =#

"""
    study_parameters(::Val{:DatabaseGenerator})

Generates a database of dds from `ini` and `act` based on ranges specified in `ini` (i.e. `ini.equilibrium.R0 = 5.0 ↔ [4.0, 10.0]`)

It's also possible to run the database generator on Vector of `ini`s and `act`s. NOTE: the length of the `ini`s and `act`s must be the same.

There is a example notebook in `FUSE_examples/study_database_generator.ipynb` that goes through the steps of setting up, running and analyzing this study
"""
function study_parameters(::Val{:DatabaseGenerator})
    return FUSEparameters__ParametersStudyDatabaseGenerator{Real}()
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
    database_policy::Switch{Symbol} = study_common_parameters(; database_policy=:single_hdf5)
end

mutable struct StudyDatabaseGenerator{T<:Real} <: AbstractStudy
    sty::OverrideParameters{T,FUSEparameters__ParametersStudyDatabaseGenerator{T}}
    ini::Union{ParametersAllInits,Vector{<:ParametersAllInits}}
    act::Union{ParametersAllActors,Vector{<:ParametersAllActors}}
    dataframe::Union{DataFrame,Missing}
    iterator::Union{Vector{Int},Missing}
    workflow::Union{Function,Missing}
end

function StudyDatabaseGenerator(sty::ParametersStudy, ini::ParametersAllInits, act::ParametersAllActors; kw...)
    sty = OverrideParameters(sty; kw...)
    study = StudyDatabaseGenerator(sty, ini, act, missing, missing, missing)
    return study
end

function StudyDatabaseGenerator(sty::ParametersStudy, inis::Vector{<:ParametersAllInits}, acts::Vector{<:ParametersAllActors}; kw...)
    @assert length(inis) == length(acts)
    sty = OverrideParameters(sty; kw...)
    if sty.n_simulations ≠ length(inis)
        @warn "sty.n_simulations is set to legth(inis)=$(length(inis))"
        sty.n_simulations = length(inis)
    end
    study = StudyDatabaseGenerator(sty, inis, acts, missing, missing, missing)

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

    @assert (sty.n_workers == 0 || sty.n_workers == length(Distributed.workers())) "The number of workers =  $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"

    if typeof(study.ini) <: ParametersAllInits && typeof(study.act) <: ParametersAllActors
        iterator = collect(1:sty.n_simulations)
    elseif typeof(study.ini) <: Vector{<:ParametersAllInits} && typeof(study.act) <: Vector{<:ParametersAllActors}
        @assert length(study.ini) == length(study.act)
        iterator = collect(1:length(study.ini))
    else
        error("DatabaseGenerator should never be here: ini and act are incompatible")
    end

    study.iterator = iterator

    # paraller run
    println("running $(sty.n_simulations) simulations with $(sty.n_workers) workers on $(sty.server)")

    if study.sty.database_policy == :separate_folders
        FUSE.ProgressMeter.@showprogress pmap(item -> run_case(study, item), iterator)
        extract_results(study)

    elseif study.sty.database_policy == :single_hdf5

        FUSE.ProgressMeter.@showprogress pmap(item -> begin
                try
                    run_case(study, item, Val(:hdf5))
                catch e
                    if isa(e, InterruptException)
                        rethrow(e)  # or handle as needed
                    end
                end
            end, iterator)

        study.dataframe = _merge_tmp_study_files(study.sty.save_folder; cleanup=true)
    else
        error("DatabaseGenerator should never be here: database_policy must be either `:separate_folders` or `:single_hdf5`")
    end

    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"
    end

    return study
end

"""
    run_case(study::AbstractStudy, item::String)

Run a single case based by setting up a dd from ini and act and then executing the workflow_DatabaseGenerator workflow (feel free to change the workflow based on your needs)
"""
function run_case(study::AbstractStudy, item::Int)
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

    # ini/act variations
    if typeof(study.ini) <: ParametersAllInits
        ini = rand(study.ini)
    elseif typeof(study.ini) <: Vector{<:ParametersAllInits}
        ini = study.ini[item]
    end
    if typeof(study.act) <: ParametersAllActors
        act = rand(study.act)
    elseif typeof(study.act) <: Vector{<:ParametersAllActors}
        act = study.act[item]
    end

    dd = IMAS.dd()

    try
        redirect_stdout(file_log)
        redirect_stderr(file_log)

        study.workflow(dd, ini, act)

        # save simulation data to directory
        save(savedir, sty.save_dd ? dd : nothing, ini, act; timer=true, freeze=false, overwrite_files=true)

        return nothing
    catch e
        if isa(e, InterruptException)
            rethrow(e)
        end

        # save empty dd and error to directory
        save(savedir, nothing, ini, act; error=e, timer=true, freeze=false, overwrite_files=true)
    finally
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        cd(original_dir)
        close(file_log)
    end
end

function run_case(study::AbstractStudy, item::Int, ::Val{:hdf5}; kw...)
    sty = study.sty
    @assert isa(study.workflow, Function) "Make sure to specicy a workflow to study.workflow that takes dd, ini , act as arguments"

    original_dir = pwd()
    if !isdir(sty.save_folder)
        mkdir(sty.save_folder)
    end
    cd(sty.save_folder)

    # Redirect stdout and stderr to the file
    original_stdout = stdout  # Save the original stdout
    original_stderr = stderr  # Save the original stderr

    # ini/act variations
    if typeof(study.ini) <: ParametersAllInits
        ini = rand(study.ini)
    elseif typeof(study.ini) <: Vector{<:ParametersAllInits}
        ini = study.ini[item]
    end

    if typeof(study.act) <: ParametersAllActors
        act = rand(study.act)
    elseif typeof(study.act) <: Vector{<:ParametersAllActors}
        act = study.act[item]
    end

    dd = IMAS.dd()


    zero_pad_length = length(string(sty.n_simulations))
    parent_group = "/case$(lpad(item, zero_pad_length, "0"))"


    tmp_log_filename = "tmp_log_worker_$(Distributed.myid())_pid$(getpid())_case_$item.txt"
    tmp_log_io = open(tmp_log_filename, "w+")

    myid = Distributed.myid()
    start_time = time()

    try
        redirect_stdout(tmp_log_io)
        redirect_stderr(tmp_log_io)

        study.workflow(dd, ini, act)

        df = DataFrame(IMAS.extract(dd, :all))
        df[!, :case] = fill(item, nrow(df))
        df[!, :dir] = fill(sty.save_folder, nrow(df))
        df[!, :gparent] = fill(parent_group, nrow(df))
        df[!, :status] = fill("success", nrow(df))
        df[!, :worker_id] = fill(myid, nrow(df))
        df[!, :elapsed_time] = fill(time() - start_time, nrow(df))

        save_study_database("tmp_h5_output", parent_group, (sty.save_dd ? dd : IMAS.dd()), ini, act, tmp_log_io; timer=true, freeze=false, overwrite_groups=true, kw...)

        # Write into temporary csv files, in case the whole Julia session is crashed
        tmp_csv_folder = "tmp_csv_output"
        if !isdir(tmp_csv_folder)
            mkdir(tmp_csv_folder)
        end
        csv_filepath = joinpath(tmp_csv_folder, "extract_success_pid$(getpid()).csv")
        if isfile(csv_filepath)
            CSV.write(csv_filepath, df; append=true, header=false)
        else
            CSV.write(csv_filepath, df)
        end

        return df

    catch e
        if isa(e, InterruptException)
            rethrow(e)
        end

        df = DataFrame(; case=item, dir=sty.save_folder, gparent=parent_group, status="fail")
        df[!, :worker_id] = fill(myid, nrow(df))
        df[!, :elapsed_time] = fill(time() - start_time, nrow(df))

        # save empty dd and error to directory
        save_study_database("tmp_h5_output", parent_group, nothing, ini, act, tmp_log_io; error_info=e, timer=true, freeze=false, overwrite_groups=true, kw...)

        # Write into temporary csv files, in case the whole Julia session is crashed
        tmp_csv_folder = "tmp_csv_output"
        if !isdir(tmp_csv_folder)
            mkdir(tmp_csv_folder)
        end
        csv_filepath = joinpath(tmp_csv_folder, "extract_fail_pid$(getpid()).csv")
        if isfile(csv_filepath)
            CSV.write(csv_filepath, df; append=true, header=false)
        else
            CSV.write(csv_filepath, df)
        end

        return df

    finally
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        close(tmp_log_io)
        rm(tmp_log_filename; force=true)

        cd(original_dir)
    end
end