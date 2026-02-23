import HDF5

#= ====================== =#
#  StudyDatabaseGenerator  #
#= ====================== =#

"""
    study_parameters(::Val{:DatabaseGenerator})

Generates a database of dds from `ini` and `act` based on ranges specified in `ini` (i.e. `ini.equilibrium.R0 = 5.0 ↔ [4.0, 10.0]`)

It's also possible to run the database generator on Vector of `ini`s and `act`s. NOTE: the length of the `ini`s and `act`s must be the same.

For large parameter sweeps (100k+ combinations), use the generator-based constructor to avoid
materializing all ini/act objects in memory. See `StudyDatabaseGenerator(sty, generator)`.

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
    ini::Union{ParametersAllInits,Vector{<:ParametersAllInits},Nothing}
    act::Union{ParametersAllActors,Vector{<:ParametersAllActors},Nothing}
    generator::Union{Nothing,Function}
    dataframe::Union{DataFrame,Missing}
    iterator::Union{Vector{Int},Missing}
    workflow::Union{Function,Missing}
end

function setup(study::StudyDatabaseGenerator)
    sty = study.sty
    check_and_create_file_save_mode(sty)
    parallel_environment(sty.server, sty.n_workers)
    return study
end

function StudyDatabaseGenerator(sty::ParametersStudy, ini::ParametersAllInits, act::ParametersAllActors; kw...)
    sty = OverrideParameters(sty; kw...)
    study = StudyDatabaseGenerator(sty, ini, act, nothing, missing, missing, missing)
    return setup(study)
end

function StudyDatabaseGenerator(sty::ParametersStudy, inis::Vector{<:ParametersAllInits}, acts::Vector{<:ParametersAllActors}; kw...)
    @assert length(inis) == length(acts)
    sty = OverrideParameters(sty; kw...)
    if sty.n_simulations ≠ length(inis)
        @warn "sty.n_simulations is set to length(inis)=$(length(inis))"
        sty.n_simulations = length(inis)
    end
    study = StudyDatabaseGenerator(sty, inis, acts, nothing, missing, missing, missing)

    check_and_create_file_save_mode(sty)

    parallel_environment(sty.server, sty.n_workers)

    return study
end

"""
    StudyDatabaseGenerator(sty::ParametersStudy, generator::Function; kw...)

Create a database generator study using a lazy generator function instead of pre-allocated
`ini`/`act` vectors. The `generator` function must accept an integer index and return a
`(ini::ParametersAllInits, act::ParametersAllActors)` tuple.

This avoids materializing all parameter combinations in memory at once, which is critical
for large parameter sweeps (100k+ combinations) where pre-allocating all `ini`/`act` objects
would exhaust available memory and cause each `pmap` worker to serialize the entire collection.

`sty.n_simulations` must be set before calling this constructor.

### Example:
```julia
grid = collect(Iterators.product(param1_values, param2_values, ...))
sty.n_simulations = length(grid)

function make_ini_act(item::Int)
    (p1, p2, ...) = grid[item]
    ini = deepcopy(ini_base)
    act = deepcopy(act_base)
    ini.some_param = p1
    ...
    return ini, act
end

study = FUSE.StudyDatabaseGenerator(sty, make_ini_act)
```
"""
function StudyDatabaseGenerator(sty::ParametersStudy, generator::Function; kw...)
    sty = OverrideParameters(sty; kw...)
    @assert sty.n_simulations > 0 "sty.n_simulations must be set before creating a generator-based study"
    study = StudyDatabaseGenerator(sty, nothing, nothing, generator, missing, missing, missing)

    check_and_create_file_save_mode(sty)

    parallel_environment(sty.server, sty.n_workers)

    return study
end

"""
    _resolve_ini_act(study::StudyDatabaseGenerator, item::Int)

Resolve the `(ini, act)` pair for a given item index, supporting generator functions,
vectors, and single ini/act with random sampling.
"""
function _resolve_ini_act(study::StudyDatabaseGenerator, item::Int)
    if study.generator !== nothing
        return study.generator(item)
    elseif study.ini isa Vector
        return study.ini[item], study.act[item]
    else
        return rand(study.ini), rand(study.act)
    end
end

"""
    _run(study::StudyDatabaseGenerator)

Runs the DatabaseGenerator with sty settings in parallel on designated cluster.

Uses lightweight pmap closures that avoid serializing the full ini/act collections to each
worker. The generator (or vector indexing) is evaluated lazily on the main process, and only
individual `(item, ini, act)` tuples are sent to workers one at a time. This means the
generator function can safely reference main-process globals without needing to be a closure.
"""
function _run(study::StudyDatabaseGenerator)
    sty = study.sty

    @assert (sty.n_workers == 0 || sty.n_workers == length(Distributed.workers())) "The number of workers =  $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"

    if study.generator !== nothing
        iterator = collect(1:sty.n_simulations)
    elseif study.ini isa ParametersAllInits && study.act isa ParametersAllActors
        iterator = collect(1:sty.n_simulations)
    elseif study.ini isa Vector{<:ParametersAllInits} && study.act isa Vector{<:ParametersAllActors}
        @assert length(study.ini) == length(study.act)
        iterator = collect(1:length(study.ini))
    else
        error("DatabaseGenerator should never be here: ini and act are incompatible")
    end

    study.iterator = iterator

    # parallel run
    println("running $(sty.n_simulations) simulations with $(sty.n_workers) workers on $(sty.server)")

    # Build a lazy generator of (item, ini, act) tuples that is iterated on the main process.
    # Only the individual tuples are serialized and sent to workers one at a time.
    # This avoids: (a) materializing all ini/act in memory, (b) serializing the full
    # study to each worker, and (c) requiring user-defined generators to be closures.
    if study.generator !== nothing
        _gen = study.generator
        work_items = ((i, _gen(i)...) for i in iterator)
    elseif study.ini isa Vector
        _inis, _acts = study.ini, study.act
        work_items = ((i, _inis[i], _acts[i]) for i in iterator)
    else
        _ini, _act = study.ini, study.act
        work_items = ((i, rand(_ini), rand(_act)) for i in iterator)
    end

    _sty = study.sty
    _wf = study.workflow

    if study.sty.database_policy == :separate_folders

        FUSE.ProgressMeter.@showprogress pmap(work_items) do (item, ini, act)
            run_case(_sty, _wf, ini, act, item)
        end
        extract_results(study)

    elseif study.sty.database_policy == :single_hdf5

        FUSE.ProgressMeter.@showprogress pmap(work_items) do (item, ini, act)
            try
                run_case(_sty, _wf, ini, act, item, Val(:hdf5))
            catch e
                isa(e, InterruptException) && rethrow(e)
                @error "run_case failed for item $item on worker $(Distributed.myid())" exception = (e, catch_backtrace())
            end
        end

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

# ================================================================ #
# run_case methods that accept (study, item) for backward compat   #
# These delegate to the standalone (sty, workflow, ini, act, item) #
# methods after resolving ini/act from the study.                  #
# ================================================================ #

"""
    run_case(study::AbstractStudy, item::Int)

Run a single case for the separate_folders policy by resolving ini/act from the study.
"""
function run_case(study::AbstractStudy, item::Int)
    sty = study.sty
    @assert isa(study.workflow, Function) "Make sure to specify a workflow to study.workflow that takes dd, ini, act as arguments"

    if study isa StudyDatabaseGenerator
        ini, act = _resolve_ini_act(study, item)
    else
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
    end

    return run_case(sty, study.workflow, ini, act, item)
end

"""
    run_case(study::AbstractStudy, item::Int, ::Val{:hdf5}; kw...)

Run a single case for the single_hdf5 policy by resolving ini/act from the study.
"""
function run_case(study::AbstractStudy, item::Int, v::Val{:hdf5}; kw...)
    sty = study.sty
    @assert isa(study.workflow, Function) "Make sure to specify a workflow to study.workflow that takes dd, ini, act as arguments"

    if study isa StudyDatabaseGenerator
        ini, act = _resolve_ini_act(study, item)
    else
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
    end

    return run_case(sty, study.workflow, ini, act, item, v; kw...)
end

# ================================================================ #
# Standalone run_case methods that accept ini/act directly.        #
# These are called from _run's pmap closures to avoid capturing    #
# the full study object (with large ini/act vectors).              #
# ================================================================ #

"""
    run_case(sty, workflow::Function, ini::ParametersAllInits, act::ParametersAllActors, item::Int)

Run a single case for the separate_folders policy with explicitly provided ini/act.
This method does not reference any large study data structures, making it safe for
use in pmap closures without serializing the entire parameter collection.
"""
function run_case(sty, workflow::Function, ini::ParametersAllInits, act::ParametersAllActors, item::Int)
    original_dir = pwd()
    savedir = abspath(joinpath(sty.save_folder, "$(item)__$(Dates.now())__$(getpid())"))
    mkdir(savedir)
    cd(savedir)

    # Redirect stdout and stderr to the file
    original_stdout = stdout  # Save the original stdout
    original_stderr = stderr  # Save the original stderr
    file_log = open("log.txt", "w")

    dd = IMAS.dd()

    try
        redirect_stdout(file_log)
        redirect_stderr(file_log)

        workflow(dd, ini, act)

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

"""
    run_case(sty, workflow::Function, ini::ParametersAllInits, act::ParametersAllActors, item::Int, ::Val{:hdf5}; kw...)

Run a single case for the single_hdf5 policy with explicitly provided ini/act.
This method does not reference any large study data structures, making it safe for
use in pmap closures without serializing the entire parameter collection.
"""
function run_case(sty, workflow::Function, ini::ParametersAllInits, act::ParametersAllActors, item::Int, ::Val{:hdf5}; kw...)
    original_dir = pwd()
    if !isdir(sty.save_folder)
        mkdir(sty.save_folder)
    end
    cd(sty.save_folder)

    # Redirect stdout and stderr to the file
    original_stdout = stdout  # Save the original stdout
    original_stderr = stderr  # Save the original stderr

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

        workflow(dd, ini, act)

        xtract = IMAS.extract(dd, :all)
        merge!(xtract, extract_ini(ini))
        merge!(xtract, extract_act(act))
        df = DataFrame(xtract)
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

        dd_empty = IMAS.dd()
        xtract = IMAS.extract(dd_empty, :all)
        merge!(xtract, extract_ini(ini))
        merge!(xtract, extract_act(act))
        df = DataFrame(xtract)
        df[!, :case] = fill(item, nrow(df))
        df[!, :dir] = fill(sty.save_folder, nrow(df))
        df[!, :gparent] = fill(parent_group, nrow(df))
        df[!, :status] = fill("fail", nrow(df))
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
