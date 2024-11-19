"""
    extract_results(simulations_path::String)

Extracts informatation from all the simulation results in simulations_path and returns them inside a dataframe (This is done safely in parallel and with checkpoints) 
Note : If you want to speed up this process incrase the number of threads
"""
function extract_results(simulations_path::String)
    function get_dataframe(file_path)
        dd = IMAS.json2imas(file_path)
        df = DataFrame(IMAS.extract(dd, :all))
        df.dir = [file_path]
        df.gen = [parse(Int,split(split(file_path,"/")[end-1],"__")[1])]
        return df
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
    if isempty(json_files)
        return DataFrames.DataFrame()
    end
    df_filler = get_dataframe(json_files[1])
    column_names = names(df_filler)
    column_types = map(col -> eltype(df_filler[!, col]), names(df_filler))
    function initialize_empty_df()
        return DataFrame([name => Vector{type}() for (name, type) in zip(column_names, column_types)])
    end

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

function extract_results(study::AbstractStudy)
    study.dataframe = extract_results(study.sty.save_folder)
end