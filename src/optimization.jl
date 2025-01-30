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
        save_dd::Bool=true)

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
    save_dd::Bool=true)

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

function _optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    x::AbstractVector,
    objective_functions::AbstractVector{<:IMAS.ObjectiveFunction},
    constraint_functions::AbstractVector{<:IMAS.ConstraintFunction},
    save_folder::AbstractString,
    generation::Int,
    save_dd::Bool=true)

    tmp = optimization_engine(ini, act, actor_or_workflow, x, objective_functions, constraint_functions, save_folder, generation, save_dd)

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
    generation_offset::Int)

    # parallel evaluation of a generation
    ProgressMeter.next!(p)
    tmp = Distributed.pmap(
        x -> _optimization_engine(ini, act, actor_or_workflow, x, objective_functions, constraint_functions, save_folder, p.counter + generation_offset, save_dd),
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
