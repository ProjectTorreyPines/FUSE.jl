using Serialization
using LaTeXStrings

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
    println("analyzing study with n_threads = $(Base.Threads.nthreads())")
    extract_optimization_results(study)
    return study
end


function extract_optimization_results(study)
    function get_dataframe(file_path)
        dd = IMAS.json2imas(file_path)
        df = DataFrame(IMAS.extract(dd, :all))
        df.dir = [file_path]
        return df
    end

    function initialize_empty_df()
        return DataFrame([name => Vector{type}() for (name, type) in zip(column_names, column_types)])
    end

    # Directory paths
    simulations_path = study.sty.save_folder
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
        processed_files_set = deserialize(checkpoint_file)
    end

    # Filter out processed files
    json_files = [file for file in all_files if !(file in processed_files_set) && endswith(file, ".json")]

    df_filler = get_dataframe(json_files[1])
    column_names = names(df_filler)
    column_types = map(col -> eltype(df_filler[!, col]), names(df_filler))
    # Number of threads
    num_threads = Base.Threads.nthreads()

    # Define batch size
    batch_size = length(json_files) > 1000 ? 1000 : Int(ceil(length(json_files) / num_threads))  # Adjust based on memory constraints
    println("analyzing $(total_files) in batches of $(batch_size) with nthreads = $(num_threads)")


    # Split files into batches
    file_batches = [json_files[i:min(i + batch_size - 1, end)] for i in 1:batch_size:length(json_files)]

    # Prepare a vector to hold the DataFrames
    data_frames = Vector{DataFrame}()

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
                serialize(pf_file, local_processed_files)
                # Store DataFrame reference
                thread_data_frames[i] = local_df
            end
        end

        # Update the global list of processed files
        for i in 1:num_threads
            # Load processed files from each thread
            pf_file = joinpath(checkpoint_dir, "pf_batch$(batch_idx)_thread$(i).jls")
            if isfile(pf_file)
                local_processed_files = deserialize(pf_file)
                union!(processed_files_set, local_processed_files)
                rm(pf_file)  # Clean up if not needed anymore
            end
        end

        # Save updated processed files set
        serialize(checkpoint_file, processed_files_set)
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
    CSV.write(joinpath(study.sty.save_folder, "output.csv"), final_df)
    return study.dataframe = final_df
end


#### Analysis


function copy_and_return_dds(df, pareto_list, laptop_dir)

    dds = IMAS.dd[]
    for par in pareto_list
        path_omega = df.dir[par]
        run(`scp -r omegae:$path_omega $latop_dir`)
        p = joinpath(latop_dir, split(df.dir[par], "/")[end])
        dd, ini, act = FUSE.load(p)
        push!(dds, dd)
    end
    return dds
end


function plot_all_traces(outputs; save=false, verbose=false)

    #    outputs = opt_run.output
    #    println("Traces for $(opt_run.name)")
    #     println("============== CONSTRAINTS ==============")

    #     for q in [(s g xtring(a),string(a)) for a in runs[1].constraint_symbols]
    #         plot_population_traces(outputs,q[1],q[2];constraint=true,save)
    #     end


    println("============== OBJECTIVES ==============")

    for q in objective_traces
        plot_population_traces(outputs, q[1], q[2]; save, verbose)
    end


    println("============== ACTUATORS ==============")


    for q in actuator_traces
        plot_population_traces(outputs, q[1], q[2]; save, verbose)
    end


    println("============== OUTPUT ==============")

    for q in other_traces
        plot_population_traces(outputs, q[1], q[2]; save, verbose)
    end

    println("============== FOR FUN ==============")


    for q in forfun_traces
        plot_population_traces(outputs, q[1], q[2]; save, verbose)
    end
end


function plot_population_traces(outputs, var_name, plot_name; constraint=false, save=false, verbose=false, cname="", xlab="Generations", ticksizes=15, labfontsize=13)
    if verbose
        @show var_name, plot_name
    end


    if var_name == "R0B0" || var_name == "B0R0"
        return
    end

    N = length(outputs[:, "Pelectric_net"])
    N_blocks = 21
    #    N_blocks = 20
    Y_points = 25
    X = LinRange(0, N, N_blocks)
    xx = [string(Int(round(i / X[end] * maximum(outputs.gen)))) for i in X]
    ar = collect(1:length(xx))
    @show length(xx)
    deleteat!(ar, collect(1:4:length(xx)))
    xx[ar] .= ""
    @show xx
    xticks = (X, xx)

    function y_auto_range(y; σ=5, N=Y_points)
        y_nonan = y[@. !isnan.(y)]
        m = Statistics.median(y_nonan)
        s = Statistics.median(Statistics.median(abs.(y_nonan .- m))) * σ
        return Y = LinRange(max(m - s, minimum(y_nonan)), min(m + s, maximum(y_nonan)), N)
    end

    if var_name == "HTS"
        y = (outputs[:, "TF_material"] .== "rebco") .+ (outputs[:, "OH_material"] .== "rebco") .* 2.0
        yname = "HTS"
        Y = LinRange(0, 3.00001, Y_points)
        p = histogram2d(y; bins=(X, Y), ylabel=yname, xlabel="Individuals")
    elseif var_name == "HTS"
        y = (outputs[:, "TF_material"] .== "rebco") .+ (outputs[:, "OH_material"] .== "rebco") .* 2.0
    elseif var_name == "thermal_cycle_type"
        y = (outputs[:, "thermal_cycle_type"] .== "rankine") .* 2.0

        yname = "thermal_cycle_type"
        Y = LinRange(0, 3.00001, Y_points)
        println(" legend ->  0 :: brayton, 2 :: rankine")
        p = histogram2d(y; bins=(X, Y), ylabel=yname, xlabel="Individuals",)

    elseif var_name == "fbs"
        y = outputs[:, "ip_bs"] ./ outputs[:, "ip"]
        yname = "bootstrap fraction"
        Y = LinRange(0, 1.0, Y_points)
        println(" legend ->  0 :: brayton, 2 :: rankine")
        p = histogram2d(y; bins=(X, Y), ylabel=yname, xlabel="Individuals",)
    elseif var_name == "gen" || var_name == "generation"
        p = false
    else
        y = outputs[:, var_name]
        yname = plot_name
        Y = y_auto_range(y)
        @show Y

        p = histogram2d(y; bins=(X, Y), ylabel=yname, xlabel="Generation", xticks=xticks,
            labelfontsize=labfontsize, xtickfontsize=ticksizes, ytickfontsize=ticksizes, colorbar_tickfontsize=ticksizes, colorbar_title="\n Population density",
            thickness_scaling=1.1, right_margin=10Plots.mm)
    end
    if constraint
        hline!([0.1]; ls=:dash, color=:cyan, lw=2.0, label="Constraint value")
    end

    if p != false
        display(p)
    end
    if p != false && save
        if var_name !== "P/R0"
            savefig(p, joinpath(picture_path, "trace_$(var_name)_$(formatted_date).pdf"))
        end

    end
    return p
end

function plot_pareto_front(
    outputs,
    x_name::String,
    y_name::String,
    color_axis_name::String;
    color_axis_name_label::String="",
    pareto_only=false,
    pareto_color=:blue,
    plt::Union{Plots.Plot,Bool}=false,
    pareto_name::String="",
    pareto_numbers=true,
    return_scatter_only=false,
    ticksizes=15,
    labfontsize=13,
    dis=true,
    kwargs...
)
    special_handeling_dict = Dict(
        "PsolR"=>(outputs[:, "Psol"] ./ outputs[:, "R0"]),
        "magnet_materials"=>((outputs[:, "TF_material"] .== "rebco") .+ (outputs[:, "OH_material"] .== "rebco")) .* 2.0,
        "B0R0" => (outputs[:, "R0"] .* outputs[:, "B0"]),
        "Psol-PLH" => (outputs[:, "Psol"] .- outputs[:, "PLH"]))

    pretty_names_dict = Dict("capital_cost"=>L"\textnormal{Capital \, cost}",
        "zeff_ped"=>L"\textbf{\textnormal{Z_{eff}}}",
        "B0R0" => L"B_{0} R_{0} [Tm]", "q95"=>L"\textbf{\textnormal{q_{95}}}")        
    x = outputs[:, x_name]
    y = outputs[:, y_name]

    if haskey(special_handeling_dict, color_axis_name)
        c = special_handeling_dict[color_axis_name]
    else
        c = outputs[:, color_axis_name]
    end

    if !isempty(color_axis_name_label)
        cname = color_axis_name_label
    elseif haskey(pretty_names_dict,color_axis_name)
        cname = pretty_names_dict[color_axis_name]
    else
        cname = color_axis_name
    end

 
    # plotting
    index = 1:length(outputs[:, x_name])
    annot = map(x -> (x, :center, 3, "courier"), index)

    xlab = haskey(pretty_names_dict, x_name) ? pretty_names_dict[x_name] : x_name
    ylab = haskey(pretty_names_dict, y_name) ? pretty_names_dict[y_name] : y_name

    P = scatter(x[index], y[index]; marker_z=c[index], colorbar_tickfontsize=labfontsize - 2,
        colorbar_titlefontsize=14, xlabel=xlab, ylabel=ylab, colorbar_title=cname,
        marker=:circle, markersize=10, markerstrokewidth=0,
        alpha=0.5, labelfontsize=labfontsize, xtickfontsize=ticksizes,
        ytickfontsize=ticksizes, right_margin=6Plots.mm,label="", kwargs...)

    if return_scatter_only
        return P
    end

    # Pareto front
    pindex = index[FUSE.pareto_front([[1 ./ x[index[k]], y[index[k]]] for k in 1:length(index)])]
    sort!(pindex; by=i -> x[i])

    if pareto_numbers
        pannot = map(x -> ("\n$x", :right, 3, "courier", :red), pindex)
        println(length(pindex))
        println(pindex)
        println(outputs.dir[pindex[1]])
    else
        pannot = nothing
    end
    if isa(plt, Plots.Plot)
        return plot!(plt, x[pindex], y[pindex]; series_annotations=pannot, color=pareto_color, label="pareto front $(pareto_name)", lw=2)
    elseif !pareto_only
        if dis
            display(plot!(P, x[pindex], y[pindex]; series_annotations=pannot, color=pareto_color, label="pareto front $(pareto_name)", lw=2, legend=:bottomright))
        else
            plot!(P, x[pindex], y[pindex]; series_annotations=pannot, color=pareto_color, label="pareto front $(pareto_name)", lw=2, legend=:bottomright)
        end
        return P
    else
        return plot(x[pindex], y[pindex]; series_annotations=pannot, color=pareto_color, label="pareto front $(pareto_name)", lw=2, xlabel=xname, ylabel=yname)
    end
end

function make_gif(quant, quant_name, quant_range)
    #    q95_selections = reverse(collect(3.5:0.5:15))
    plots = []
    for q in quant_range
        println("$q")
        outputs_filtered2 = outputs_filtered[filter(i -> outputs_filtered[i, quant] < q, 1:length(outputs_filtered.gen)), :]
        #        @show outputs_filtered2.gen
        push!(plots, plot_pareto_front(outputs_filtered2, "gen", "Generation"; xxlim=(2.0, 10.0), yylim=(0.0, 22.0), pareto_numbers=false, clim=(1, 100), ticksizes=13);)
    end

    return plots
end




pretty_names_dict = Dict("capital_cost"=>L"\textnormal{Capital \, cost}",
    "zeff_ped"=>L"\textbf{\textnormal{Z_{eff}}}",
    "B0R0" => L"B_{0} R_{0} [Tm]", "q95"=>L"\textbf{\textnormal{q_{95}}}",
    
    )