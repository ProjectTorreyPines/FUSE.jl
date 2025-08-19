using Serialization: Serialization

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
    restart_workers_after_n_generations::Entry{Int} =
        Entry{Int}("-", "Runs the optimization in a safe way restarting the workers every N generations, default is never restart"; default=0)
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the database runs into")
    # Optimization related parameters
    population_size::Entry{Int} = Entry{Int}("-", "Number of individuals in a generation")
    number_of_generations::Entry{Int} = Entry{Int}("-", "Number generations")
    database_policy::Switch{Symbol} = study_common_parameters(; database_policy=:single_hdf5)
end

mutable struct StudyMultiObjectiveOptimizer{T<:Real} <: AbstractStudy
    sty::OverrideParameters{T,FUSEparameters__ParametersStudyMultiObjectiveOptimizer{T}}
    ini::ParametersAllInits
    act::ParametersAllActors
    constraint_functions::Vector{IMAS.ConstraintFunction}
    objective_functions::Vector{IMAS.ObjectiveFunction}
    state::Union{Nothing,Metaheuristics.State}
    dataframe::Union{DataFrame,Missing}
    datafame_filtered::Union{DataFrame,Missing}
    generation::Int
    workflow::Union{Function,Missing}
end

function StudyMultiObjectiveOptimizer(
    sty::ParametersStudy,
    ini::ParametersAllInits,
    act::ParametersAllActors,
    constraint_functions::Vector{IMAS.ConstraintFunction},
    objective_functions::Vector{IMAS.ObjectiveFunction};
    kw...
)
    sty = OverrideParameters(sty; kw...)
    study = StudyMultiObjectiveOptimizer(sty, ini, act, constraint_functions, objective_functions, nothing, missing, missing, 0, missing)
    return setup(study)
end

function _setup(study::StudyMultiObjectiveOptimizer)
    sty = study.sty

    check_and_create_file_save_mode(sty)

    parallel_environment(sty.server, sty.n_workers)

    # import FUSE and IJulia on workers
    if isdefined(Main, :IJulia)
        code = """
        using Distributed
        @everywhere import FUSE
        @everywhere import IJulia
        """
    else
        code = """
        using Distributed
        @everywhere import FUSE
        """
    end
    Base.include_string(Main, code)
    return study
end

function _run(study::StudyMultiObjectiveOptimizer)
    sty = study.sty

    @assert sty.n_workers == length(Distributed.workers()) "The number of workers =  $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"
    @assert iseven(sty.population_size) "Population size must be even"

    if sty.restart_workers_after_n_generations > 0
        # if restart_workers_after_n_generations we are going to call _run again with modified
        max_gens_per_iteration = sty.restart_workers_after_n_generations
        steps = Int(ceil(sty.number_of_generations / max_gens_per_iteration))
        sty_bkp = deepcopy(sty)
        for i in 1:steps
            try
                println("Running $max_gens_per_iteration generations ($i / $steps)")
                gens = max_gens_per_iteration
                if i == steps && mod(sty.number_of_generations, max_gens_per_iteration) != 0
                    gens = mod(sty.restart_workers_after_n_generations, max_gens_per_iteration)
                end
                sty.restart_workers_after_n_generations = 0
                sty.release_workers_after_run = false
                sty.file_save_mode = :append
                sty.number_of_generations = gens

                if study.state !== nothing
                    study.state = load_optimization(joinpath(sty.save_folder, "results.jls")).state
                    study.generation += study.state.iteration
                end
                run(study)
            catch e
                if isa(e, InterruptException)
                    rethrow(e)
                end
                @warn "error occured in step $i \n $(string(e))"
            finally
                sty.restart_workers_after_n_generations = sty_bkp.restart_workers_after_n_generations
                sty.release_workers_after_run = sty_bkp.release_workers_after_run
                sty.file_save_mode = sty_bkp.file_save_mode
            end
        end
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"

    else
        setup(study)
        optimization_parameters = Dict(
            :N => sty.population_size,
            :iterations => sty.number_of_generations,
            :continue_state => study.state,
            :save_folder => sty.save_folder)

        @assert !isempty(sty.save_folder) "Specify where you would like to store your optimization results in sty.save_folder"

        if ismissing(study.workflow)
            study.workflow = optimization_workflow_default
        end
        study.state = workflow_multiobjective_optimization(
            study.ini, study.act, study.workflow, study.objective_functions, study.constraint_functions;
            optimization_parameters..., generation_offset=study.generation, sty.database_policy,
            sty.number_of_generations, sty.population_size)

        save_optimization(
            joinpath(sty.save_folder, "optimization_state.bson"),
            study.state,
            study.ini,
            study.act,
            study.objective_functions,
            study.constraint_functions)

        if study.sty.database_policy == :separate_folders
            analyze(study; extract_results=true)
        else
            study.dataframe = _merge_tmp_study_files(sty.save_folder; cleanup=true)
            analyze(study; extract_results=false)
        end

        # Release workers after run
        if sty.release_workers_after_run
            Distributed.rmprocs(Distributed.workers())
            @info "released workers"
        end

        return study
    end
end

function _merge_tmp_study_files(save_folder::AbstractString; cleanup::Bool=false)
    @assert isdir(save_folder) "The folder (\"$save_folder\") you are trying to merge does not exist"

    merged_hdf5_filepath = joinpath(save_folder, "database.h5")

    @info "Merging temporary study files into \"$(merged_hdf5_filepath)\"..."

    # read csv files
    tmp_csv_folder = joinpath(save_folder, "tmp_csv_output")
    if !isdir(tmp_csv_folder)
        @warn "Nothing to merge in \"$save_folder\" (pwd=$(pwd()))"

        df = DataFrame()
        if isfile(merged_hdf5_filepath)
            HDF5.h5open(merged_hdf5_filepath, "r+") do fid
                if haskey(fid, "/extract.csv")
                    df = coalesce.(CSV.read(IOBuffer(fid["/extract.csv"][]), DataFrame), NaN)
                end
            end
        end
        return df
    end

    csv_files = readdir(tmp_csv_folder; join=true)
    dfs = [CSV.read(file, DataFrame) for file in csv_files]

    merged_df = reduce(vcat, dfs; cols=:union)
    sort!(merged_df, "gparent")

    if "gen" in names(merged_df)
        # This is a multi-objective optimization
        leading_cols = ["gparent", "status", "gen", "case", "Ngen", "Ncase", "dir", "worker_id", "elapsed_time"]
        h5_group_search_depth = 2
    else
        # This is a database generator
        leading_cols = ["gparent", "status", "case", "dir", "worker_id", "elapsed_time"]
        h5_group_search_depth = 1
    end
    remaining = setdiff(names(merged_df), leading_cols)
    desired_order = vcat(leading_cols, sort(remaining))
    merged_df = coalesce.(merged_df[:, desired_order], NaN)

    IMAS.h5merge(merged_hdf5_filepath,
        joinpath(save_folder, "tmp_h5_output");
        h5_group_search_depth=h5_group_search_depth,
        h5_strip_group_prefix=true,
        cleanup)

    # Add merged_df into the mergedh5 file
    HDF5.h5open(merged_hdf5_filepath, "r+") do fid
        # Check if extract.csv already exists and append to it
        if haskey(fid, "/extract.csv")
            existing_df = coalesce.(CSV.read(IOBuffer(fid["/extract.csv"][]), DataFrame), NaN)
            merged_df = vcat(existing_df, merged_df; cols=:union)

            HDF5.delete_object(fid, "/extract.csv")
        end

        unique!(merged_df)
        sort!(merged_df, "gparent")

        io_buffer = IOBuffer()
        CSV.write(io_buffer, merged_df)
        csv_text = String(take!(io_buffer))
        HDF5.write(fid, "extract.csv", csv_text)
        attr = HDF5.attrs(fid["/extract.csv"])
        attr["date_time"] = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        return attr["FUSE_version"] = string(pkgversion(FUSE))
    end

    merged_csv_filepath = joinpath(save_folder, "extract.csv")
    CSV.write(merged_csv_filepath, merged_df)

    if cleanup
        rm(tmp_csv_folder; recursive=true)
    end

    return merged_df
end

"""
    _analyze(study::StudyMultiObjectiveOptimizer; extract_results::Bool=true)
"""
function _analyze(study::StudyMultiObjectiveOptimizer; extract_results::Bool=true)
    if extract_results
        extract_results(study)
    end
    if !isempty(study.dataframe)
        study.datafame_filtered = filter_outputs(study.dataframe, [o.name for o in study.objective_functions])
    end
    return study
end

"""
    filter_outputs(outputs::DataFrame,constraint_symbols::Vector{Symbol})

Filters the dataframe to the constraints you pass.

Common usage will be `df_filtered = FUSE.filter_outputs(df, constraint_list)`
"""
function filter_outputs(outputs::DataFrame, constraint_symbols::Vector{Symbol})
    n = nrow(outputs)
    constraint_values = [outputs[i, key] for key in constraint_symbols, i in 1:n]
    all_constraint_idxs = findall(i -> all(x -> x == 0.0, constraint_values[:, i]), 1:n)
    return outputs[all_constraint_idxs, :]
end

"""
    optimization_workflow_default(ini::ParametersAllInits, act::ParametersAllActors)

Default optimization workflow when study.workflow isn't set, initializes and runs the whole facility actor
"""
function optimization_workflow_default(ini::ParametersAllInits, act::ParametersAllActors)
    dd = FUSE.IMAS.dd()
    FUSE.init(dd, ini, act)
    FUSE.ActorFluxMatcher(dd, act)
    return dd
end