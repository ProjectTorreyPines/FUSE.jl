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
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the database runs into")
    # Optimization related parameters
    population_size::Entry{Int} = Entry{Int}("-", "Number of individuals in a generation")
    number_of_generations::Entry{Int} = Entry{Int}("-", "Number generations")
end

mutable struct StudyMultiObjectiveOptimizer <: AbstractStudy
    sty::FUSEparameters__ParametersStudyMultiObjectiveOptimizer
    ini::ParametersAllInits
    act::ParametersAllActors
    constraint_functions::Vector{IMAS.ConstraintFunction}
    objective_functions::Vector{IMAS.ObjectiveFunction}
    dataframe::Union{DataFrame,Missing}
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
    @assert !isempty(sty.save_folder) "Specify where you would like to store your optimization results in sty.save_folder"

    optimization_parameters = Dict(
        :N => sty.population_size,
        :iterations => sty.number_of_generations,
        :continue_state => nothing,
        :save_folder => sty.save_folder)
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

    analyze(study)
    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"
    end

    return study
end

function _analyze(study::StudyMultiObjectiveOptimizer)
    println("Analyzing study with n_threads = $(Base.Threads.nthreads())")
    extract_optimization_results(study)
    return study
end

"""
    extract_optimization_results(simulations_path::String)

Extracts informatation from all the simulation results in simulations_path and returns them inside a dataframe (This is done safely in parallel and with checkpoints) 
Note : If you want to speed up this process incrase the number of threads
"""
function extract_optimization_results(simulations_path::String)
    function get_dataframe(file_path)
        dd = IMAS.json2imas(file_path)
        df = DataFrame(IMAS.extract(dd, :all))
        df.dir = [file_path]
        df.gen = [parse(Int,split(file_path,"/")[end-1])]
        return df
    end

    df_filler = get_dataframe(json_files[1])
    column_names = names(df_filler)
    column_types = map(col -> eltype(df_filler[!, col]), names(df_filler))
    function initialize_empty_df()
        return DataFrame([name => Vector{type}() for (name, type) in zip(column_names, column_types)])
    end

    # Directory paths
    checkpoint_dir = joinpath(simulations_path, "checkpoints")

    # Ensure the checkpoint directory exists
    if !isdir(checkpoint_dir)
        mkpath(checkpoint_dir)
    end

    # Get list of JSON files
    all_files = filter(isdir, sort!(readdir(simulations_path; join=true)))

    dirs = sort(filter(x -> !isfile(joinpath(x, "error.txt")) && isfile(joinpath(x, "dd.json")), all_files); by=x -> parse(Int, (split(split(x, "/")[end], "__")[1])))
    all_files = [joinpath(dir, "dd.json") for dir in dirs]
    total_files = length(all_files)

    # Load list of already processed files if exists
    processed_files_set = Set{String}()
    checkpoint_file = joinpath(checkpoint_dir, "processed_files.jls")
    if isfile(checkpoint_file)
        processed_files_set = Serialization.deserialize(checkpoint_file)
    end

    # Filter out processed files
    json_files = [file for file in all_files if !(file in processed_files_set) && endswith(file, ".json")]

    # Number of threads
    num_threads = Base.Threads.nthreads()

    # Define batch size
    batch_size = length(json_files) > 1000 ? 1000 : Int(ceil(length(json_files) / num_threads))  # Adjust based on memory constraints
    println("analyzing $(total_files) in batches of $(batch_size) with nthreads = $(num_threads)")

    # Split files into batches
    file_batches = [json_files[i:min(i + batch_size - 1, end)] for i in 1:batch_size:length(json_files)]

    # Process batches
    for (batch_idx, batch_files) in enumerate(file_batches)
        println("Processing batch $(batch_idx) of $(length(file_batches))")

        # Split batch files among threads
        file_chunks = [batch_files[i:num_threads:end] for i in 1:num_threads]

        # Prepare a vector to hold DataFrames from threads
        thread_data_frames = Vector{DataFrame}(undef, num_threads)

        # 
        @sync for i in 1:num_threads
            @async begin
                local_df = initialize_empty_df()
                local_processed_files = String[]
                for file in file_chunks[i]
                    try
                        df = get_dataframe(file)
                        append!(local_df, df)
                        push!(local_processed_files, file)
                    catch e
                        if isa(e, InterruptException)
                            rethrow(e)
                        end
                        @warn "Failed to process $file: $e"
                    end
                end
                # Save local DataFrame to disk
                df_file = joinpath(checkpoint_dir, "df_batch$(batch_idx)_thread$(i).csv")
                CSV.write(df_file, local_df)
                # Save list of processed files
                pf_file = joinpath(checkpoint_dir, "pf_batch$(batch_idx)_thread$(i).jls")
                Serialization.serialize(pf_file, local_processed_files)
                # Store DataFrame reference
                thread_data_frames[i] = local_df
            end
        end

        # Update the global list of processed files
        for i in 1:num_threads
            # Load processed files from each thread
            pf_file = joinpath(checkpoint_dir, "pf_batch$(batch_idx)_thread$(i).jls")
            if isfile(pf_file)
                local_processed_files = Serialization.deserialize(pf_file)
                union!(processed_files_set, local_processed_files)
                rm(pf_file)  # Clean up if not needed anymore
            end
        end

        # Save updated processed files set
        Serialization.serialize(checkpoint_file, processed_files_set)
    end

    final_df = DataFrame()
    for batch_idx in 1:length(file_batches)
        for i in 1:num_threads
            df_file = joinpath(checkpoint_dir, "df_batch$(batch_idx)_thread$(i).csv")
            if isfile(df_file)
                df = CSV.read(df_file, DataFrame)
                append!(final_df, df)
                rm(df_file)  # Clean up if not needed anymore
            end
        end
    end

    # Export the final DataFrame
    CSV.write(joinpath(simulations_path, "output.csv"), final_df)
    return final_df
end

function extract_optimization_results(study::StudyMultiObjectiveOptimizer)
    study.dataframe = extract_optimization_results(study.sty.save_folder)
end

"""
    filter_outputs(outputs::DataFrame,constraint_symbols::Vector{Symbol})

Filters the dataframe to the constraints you pass.
Common usage will be df_filtered = FUSE.filter_outputs(df, constraint_list)
"""
function filter_outputs(outputs::DataFrame,constraint_symbols::Vector{Symbol})
    n = length(outputs.Pelectric_net)
    constraint_values = [ outputs[i,key] for key in constraint_symbols, i in 1:n]
    all_constraint_idxs = findall(i -> all(x -> x == 0.0, constraint_values[:,i]),1:n)
    return outputs[all_constraint_idxs,:]
end