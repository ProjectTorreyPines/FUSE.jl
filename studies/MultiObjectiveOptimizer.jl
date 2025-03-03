import Serialization
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
    database_policy::Switch{Symbol} = study_common_parameters(; database_policy=:separate_folders)
end

mutable struct StudyMultiObjectiveOptimizer <: AbstractStudy
    sty::FUSEparameters__ParametersStudyMultiObjectiveOptimizer
    ini::ParametersAllInits
    act::ParametersAllActors
    constraint_functions::Vector{IMAS.ConstraintFunction}
    objective_functions::Vector{IMAS.ObjectiveFunction}
    state::Union{Nothing,Metaheuristics.State}
    dataframe::Union{DataFrame,Missing}
    datafame_filtered::Union{DataFrame,Missing}
    generation::Int
end

function StudyMultiObjectiveOptimizer(
    sty::ParametersStudy,
    ini::ParametersAllInits,
    act::ParametersAllActors,
    constraint_functions::Vector{IMAS.ConstraintFunction},
    objective_functions::Vector{IMAS.ObjectiveFunction};
    kw...
)
    sty = sty(kw...)
    study = StudyMultiObjectiveOptimizer(sty, ini, act, constraint_functions, objective_functions, nothing, missing, missing, 0)
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
                    gens = mod(gen, max_gens_per_iteration)
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
        state = workflow_multiobjective_optimization(
            study.ini, study.act, ActorWholeFacility, study.objective_functions, study.constraint_functions;
            optimization_parameters..., generation_offset=study.generation, database_policy=sty.database_policy,
            number_of_generations = sty.number_of_generations, population_size = sty.population_size)

        study.state = state

        save_optimization(
            joinpath(sty.save_folder, "optimization_state.bson"),
            state,
            study.ini,
            study.act,
            study.objective_functions,
            study.constraint_functions)

        if study.sty.database_policy == :separate_folders
            analyze(study; extract_results=true)
        else
            _merge_tmp_files(study; cleanup=true)
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

function _merge_tmp_files(study::StudyMultiObjectiveOptimizer; cleanup::Bool=false)
    date_time_str = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    merged_hdf5_filename = "generations_$(date_time_str).h5"

    save_folder = study.sty.save_folder

    IMAS.h5merge(joinpath(save_folder, merged_hdf5_filename),
        joinpath(save_folder, "tmp_h5_output");
        h5_group_search_depth=2,
        h5_strip_group_prefix=true,
        cleanup)

    # read csv files
    tmp_csv_folder = joinpath(save_folder, "tmp_csv_output")
    csv_files = readdir(tmp_csv_folder; join=true)
    dfs = [CSV.read(file, DataFrame) for file in csv_files]

    study.dataframe = reduce(vcat, dfs)

    merged_csv_filepath = joinpath(save_folder, "extract_$(date_time_str).csv")
    CSV.write(merged_csv_filepath, study.dataframe)

    # write study.dataframe into the mergedh5 file as well
    HDF5.h5open(joinpath(save_folder, merged_hdf5_filename), "r+") do fid
        io_buffer = IOBuffer()
        CSV.write(io_buffer, study.dataframe)
        csv_text = String(take!(io_buffer))
        HDF5.write(fid, "extract.csv", csv_text)
        attr = HDF5.attrs(fid["/extract.csv"])
        attr["date_time"] = date_time_str
        return attr["FUSE_version"] = string(pkgversion(FUSE))
    end

    if cleanup
        rm(tmp_csv_folder; recursive=true)
    end
end


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
Common usage will be df_filtered = FUSE.filter_outputs(df, constraint_list)
"""
function filter_outputs(outputs::DataFrame, constraint_symbols::Vector{Symbol})
    n = length(outputs.Pelectric_net)
    constraint_values = [outputs[i, key] for key in constraint_symbols, i in 1:n]
    all_constraint_idxs = findall(i -> all(x -> x == 0.0, constraint_values[:, i]), 1:n)
    return outputs[all_constraint_idxs, :]
end