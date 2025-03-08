import Metaheuristics
import Distributed
import Dates

# =================== #
# Optimization engine #
# =================== #
"""
    optimization_engine(
        ini::ParametersAllInits,
        act::ParametersAllActors,
        actor_or_workflow::Union{Type{<:AbstractActor},Function},
        x::AbstractVector,
        objective_functions::AbstractVector{<:IMAS.ObjectiveFunction},
        constraint_functions::AbstractVector{<:IMAS.ConstraintFunction},
        save_folder::AbstractString,
        generation::Int,
        save_dd::Bool)

NOTE: This function is run by the worker nodes
"""
function optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    x::AbstractVector,
    objective_functions::AbstractVector{<:IMAS.ObjectiveFunction},
    constraint_functions::AbstractVector{<:IMAS.ConstraintFunction},
    save_folder::AbstractString,
    generation::Int,
    save_dd::Bool)

    # create working directory
    original_dir = pwd()
    savedir = abspath(joinpath(save_folder, "$(generation)__$(join(split(string(Dates.now()),":"),"-"))__$(getpid())"))
    mkdir(savedir)
    cd(savedir)

    # Redirect stdout and stderr to the file
    original_stdout = stdout  # Save the original stdout
    original_stderr = stderr  # Save the original stderr
    file_log = open("log.txt", "w")

    try
        redirect_stdout(file_log)
        redirect_stderr(file_log)

        # deepcopy ini/act to avoid changes
        ini = deepcopy(ini)
        act = deepcopy(act)

        # update ini based on input optimization vector `x`
        @show x
        parameters_from_opt!(ini, x)
        @show x

        # attempt to release memory
        malloc_trim_if_glibc()

        # run the problem
        if typeof(actor_or_workflow) <: Function
            dd = actor_or_workflow(ini, act)
        else
            dd = init(ini, act)
            actor = actor_or_workflow(dd, act)
            dd = actor.dd
        end

        # save simulation data to directory
        save(savedir, save_dd ? dd : nothing, ini, act; timer=true, freeze=false, overwrite_files=true)

        # evaluate multiple objectives
        ff = collect(map(f -> nan2inf(f(dd)), objective_functions))
        gg = collect(map(g -> nan2inf(g(dd)), constraint_functions))
        hh = Float64[]

        println("finished")
        @show ff
        @show gg
        @show hh

        return ff, gg, hh

    catch error
        if isa(error, InterruptException)
            rethrow(error)
        end

        # save empty dd and error to directory
        save(savedir, nothing, ini, act; error, timer=true, freeze=false, overwrite_files=true)

        # rethrow(e) # uncomment for debugging purposes

        ff = Float64[Inf for f in objective_functions]
        gg = Float64[Inf for g in constraint_functions]
        hh = Float64[]

        println("failed")
        @show ff
        @show gg
        @show hh

        return ff, gg, hh

    finally
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        cd(original_dir)
        close(file_log)
    end
end

"""
    optimization_engine(
        ini::ParametersAllInits,
        act::ParametersAllActors,
        actor_or_workflow::Union{Type{<:AbstractActor},Function},
        x::AbstractVector,
        objective_functions::AbstractVector{<:IMAS.ObjectiveFunction},
        constraint_functions::AbstractVector{<:IMAS.ConstraintFunction},
        save_folder::AbstractString,
        generation::Int,
        save_dd::Bool,
        ::Type{Val{:hdf5}};
        case_index::Union{Nothing,Int}=nothing
        )

NOTE: This function is run by the worker nodes
"""
function optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    x::AbstractVector,
    objective_functions::AbstractVector{<:IMAS.ObjectiveFunction},
    constraint_functions::AbstractVector{<:IMAS.ConstraintFunction},
    save_folder::AbstractString,
    generation::Int,
    save_dd::Bool,
    ::Type{Val{:hdf5}};
    case_index::Union{Nothing,Int}=nothing,
    kw...
)

    # create working directory
    original_dir = pwd()
    if !isdir(save_folder)
        mkdir(save_folder)
    end
    cd(save_folder)

    file_lock_channels = get(kw, :file_lock_channels, nothing)

    number_of_generations = get(kw, :number_of_generations, 10000)
    population_size = get(kw, :population_size, 10000)

    # Redirect stdout and stderr to the file
    original_stdout = stdout  # Save the original stdout
    original_stderr = stderr  # Save the original stderr

    if isnothing(case_index)
        parent_group = "/gen$(lpad(generation,Lpad_gen,"0"))/pid$(getpid())"
        tmp_log_filename = "tmp_log_pid$(getpid()).txt"
        tmp_log_io = open(joinpath(tmp_log_folder, "pid$(getpid()).txt"), "w+")
    else
        Lpad_gen = length(string(number_of_generations))
        Lpad_case = length(string(population_size))
        parent_group = "/gen$(lpad(generation,Lpad_gen,"0"))/case$(lpad(case_index,Lpad_case, "0"))"
        tmp_log_filename = "tmp_log_worker_$(Distributed.myid())_pid$(getpid())_gen_$(generation)_case_$case_index.txt"
    end
    tmp_log_io = open(tmp_log_filename, "w+")

    myid = Distributed.myid()

    try
        redirect_stdout(tmp_log_io)
        redirect_stderr(tmp_log_io)

        # deepcopy ini/act to avoid changes
        ini = deepcopy(ini)
        act = deepcopy(act)        # update ini based on input optimization vector `x`
        @show x
        parameters_from_opt!(ini, x)
        @show x

        # attempt to release memory
        malloc_trim_if_glibc()

        # run the problem
        if typeof(actor_or_workflow) <: Function
            dd = actor_or_workflow(ini, act)
        else
            dd = init(ini, act)
            actor = actor_or_workflow(dd, act)
            dd = actor.dd
        end

        wait_for_unlock!(file_lock_channels.m2w[myid])
        put!(file_lock_channels.w2m[myid], :lock)

        # save simulation data
        save_database("tmp_h5_output", parent_group, (save_dd ? dd : nothing), ini, act, tmp_log_io;
            timer=true, freeze=false, overwrite_groups=true)

        df = DataFrame(IMAS.extract(dd, :all))
        df[!, :dir] = [relpath(".", original_dir)]
        df[!, :gen] = fill(generation, nrow(df))
        df[!, :case] = fill(case_index, nrow(df))
        df[!, :gparent] = fill(parent_group, nrow(df))
        df[!, :Ngen] = fill(number_of_generations, nrow(df))
        df[!, :Ncase] = fill(population_size, nrow(df))
        df[!, :status] = fill("success", nrow(df))

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
        sleep(1) # wait a bit to make sure the file is written
        put!(file_lock_channels.w2m[myid], :unlock)

        # evaluate multiple objectives
        ff = collect(map(f -> nan2inf(f(dd)), objective_functions))
        gg = collect(map(g -> nan2inf(g(dd)), constraint_functions))
        hh = Float64[]

        println("finished")
        @show ff
        @show gg
        @show hh

        return ff, gg, hh

    catch error
        if isa(error, InterruptException)
            rethrow(error)
        end
        wait_for_unlock!(file_lock_channels.m2w[myid])
        put!(file_lock_channels.w2m[myid], :lock)

        # save empty dd and error to directory
        save_database("tmp_h5_output", parent_group, nothing, ini, act, tmp_log_io;
            error_info=error, timer=true, freeze=false, overwrite_groups=true, kw...)

        df = DataFrame()
        df[!, :dir] = [relpath(".", original_dir)]
        df[!, :gen] = fill(generation, nrow(df))
        df[!, :case] = fill(case_index, nrow(df))
        df[!, :gparent] = fill(parent_group, nrow(df))
        df[!, :Ngen] = fill(number_of_generations, nrow(df))
        df[!, :Ncase] = fill(population_size, nrow(df))
        df[!, :status] = fill("fail", nrow(df))

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
        sleep(1) # wait a bit to make sure the file is written
        put!(file_lock_channels.w2m[myid], :unlock)

        # rethrow(e) # uncomment for debugging purposes

        ff = Float64[Inf for f in objective_functions]
        gg = Float64[Inf for g in constraint_functions]
        hh = Float64[]

        println("failed")
        @show ff
        @show gg
        @show hh

        return ff, gg, hh

    finally

        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        close(tmp_log_io)
        rm(tmp_log_filename; force=true)

        cd(original_dir)
    end
end


function _optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    x::AbstractVector,
    objective_functions::AbstractVector{<:IMAS.ObjectiveFunction},
    constraint_functions::AbstractVector{<:IMAS.ConstraintFunction},
    save_folder::AbstractString,
    generation::Int,
    save_dd::Bool;
    case_index::Union{Nothing, Int}=nothing,
    kw...)

    database_policy = get(kw, :database_policy, :separate_folders)

    if database_policy == :separate_folders
        tmp = optimization_engine(ini, act, actor_or_workflow, x, objective_functions, constraint_functions, save_folder, generation, save_dd)
    elseif database_policy == :single_hdf5
        tmp = optimization_engine(ini, act, actor_or_workflow, x, objective_functions, constraint_functions, save_folder, generation, save_dd, Val{:hdf5}; case_index, kw...)
    else
        error("database_policy must be `separate_folders` or `:single_hdf5`")
    end

    GC.gc()
    return tmp
end

"""
    optimization_engine(
        ini::ParametersAllInits,
        act::ParametersAllActors,
        actor_or_workflow::Union{Type{<:AbstractActor},Function},
        X::AbstractMatrix,
        objective_functions::AbstractVector{<:IMAS.ObjectiveFunction},
        constraint_functions::AbstractVector{<:IMAS.ConstraintFunction},
        save_folder::AbstractString,
        save_dd::Bool,
        p::ProgressMeter.Progress)

NOTE: this function is run by the master process
"""
function optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    X::AbstractMatrix,
    objective_functions::AbstractVector{<:IMAS.ObjectiveFunction},
    constraint_functions::AbstractVector{<:IMAS.ConstraintFunction},
    save_folder::AbstractString,
    save_dd::Bool,
    p::ProgressMeter.Progress,
    generation_offset::Int;
    kw...)

    # parallel evaluation of a generation
    ProgressMeter.next!(p)
    tmp = Distributed.pmap(
        (k,x) -> _optimization_engine(ini, act, actor_or_workflow, x, objective_functions, constraint_functions, save_folder, p.counter + generation_offset, save_dd; case_index=k, kw...),
        1:size(X,1),
        [X[k, :] for k in 1:size(X)[1]]
    )
    F = zeros(size(X)[1], length(tmp[1][1]))
    G = zeros(size(X)[1], max(length(tmp[1][2]), 1))
    H = zeros(size(X)[1], max(length(tmp[1][3]), 1))
    for k in 1:size(X)[1]
        f, g, h = tmp[k]
        F[k, :] .= f
        if !isempty(g)
            G[k, :] .= g
        end
        if !isempty(h)
            H[k, :] .= h
        end
    end
    return F, G, H
end

"""
    nan2inf(x::Float64)

Turn NaNs into Inf
"""
function nan2inf(x::Float64)
    if isnan(x)
        return Inf
    else
        return x
    end
end
