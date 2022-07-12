import Distributed
import BSON

mutable struct MultiobjectiveOptimizationResults
    workflow::Union{DataType,Function}
    ini::ParametersInit
    act::ParametersActor
    state::Metaheuristics.State
    opt_ini::Vector{<:Parameter}
    objectives_functions::Vector{<:ObjectiveFunction}
end

"""
    workflow_multiobjective_optimization(
        ini::ParametersInit,
        act::ParametersActor,
        actor_or_workflow::Union{DataType, Function},
        objectives_functions::Vector{<:ObjectiveFunction}=ObjectiveFunction[];
        N::Int=10,
        iterations::Int=N,
        continue_results::Union{Missing,MultiobjectiveOptimizationResults}=missing)

Multi-objective optimization of either an `actor(dd, act)` or a `workflow(ini, act)`
"""
function workflow_multiobjective_optimization(
    ini::ParametersInit,
    act::ParametersActor,
    actor_or_workflow::Union{DataType,Function},
    objectives_functions::Vector{<:ObjectiveFunction}=ObjectiveFunction[];
    N::Int=10,
    iterations::Int=N,
    continue_results::Union{Missing,MultiobjectiveOptimizationResults}=missing
)

    if mod(N, 2) > 0
        error("workflow_multiobjective_optimization population size `N` must be an even number")
    end

    println("Running on $(nprocs()-1) worker processes")
    if isempty(objectives_functions)
        error(
            "Must specify objective functions. Available pre-baked functions from ObjectivesFunctionsLibrary:\n  * " *
            join(keys(ObjectivesFunctionsLibrary), "\n  * "),
        )
    end

    # itentify optimization variables in ini
    opt_ini = opt_parameters(ini)
    println("== Actuators ==")
    for optpar in opt_ini
        println(optpar)
    end
    println()
    println("== Objectives ==")
    for objf in objectives_functions
        println(objf)
    end
    println()

    # optimization boundaries
    bounds = [[optpar.lower for optpar in opt_ini] [optpar.upper for optpar in opt_ini]]'

    # test running function once with nominal parameters useful to catch bugs quickly.
    # Use @everywhere to trigger compilation on all worker nodes.
    if typeof(actor_or_workflow) <: DataType
        actor_or_workflow(init(ini, act), act)
    else
        actor_or_workflow(ini, act)
    end

    # optimize
    options = Metaheuristics.Options(seed=1, parallel_evaluation=true, store_convergence=true, iterations=iterations)
    algorithm = Metaheuristics.NSGA2(; N, options)
    if continue_results !== missing
        println("Restarting simulation")
        algorithm.status = continue_results.state
    end
    flush(stdout)
    p = Progress(iterations; desc="Iteration", showspeed=true)
    @time state = Metaheuristics.optimize(X -> optimization_engine(ini, act, actor_or_workflow, X, opt_ini, objectives_functions, p), bounds, algorithm)

    return MultiobjectiveOptimizationResults(actor_or_workflow, ini, act, state, opt_ini, objectives_functions)
end

"""
    save_optimization(filename::AbstractString, results::FUSE.MultiobjectiveOptimizationResults)

Save MultiobjectiveOptimizationResults to file in BSON format
"""
function save_optimization(filename::AbstractString, results::FUSE.MultiobjectiveOptimizationResults)
    return BSON.bson(filename, Dict("results" => results))
end

"""
    load_optimization(filename::AbstractString)

Load MultiobjectiveOptimizationResults from file
"""
function load_optimization(filename::AbstractString)
    return BSON.load(filename, FUSE)["results"]
end

function pretty_label(objective_function::ObjectiveFunction, units="")
    txt = join(split(string(objective_function.name), "_")[2:end], " ")
    if length(units) > 0
        txt *= " [$units]"
    elseif length(objective_function.units) > 0
        txt *= " [$(objective_function.units)]"
    end
    return txt
end

function pretty_label(parameter::Parameter, units="")
    txt = replace(string(parameter._name), "_" => " ")
    if length(units) > 0
        txt *= " [$units]"
    elseif length(parameter.units) > 0
        txt *= " [$(parameter.units)]"
    end
    return txt
end

@recipe function plot_MultiobjectiveOptimizationResults(
    results::MultiobjectiveOptimizationResults,
    indexes::AbstractVector{<:Integer}=[1, 2, 3];
    color_by=0,
    design_space=false,
    pareto=true,
    max_samples=nothing,
    iterations=nothing)

    @assert length(indexes) <= 3 "plot_MultiobjectiveOptimizationResults: Cannot visualize more than 3 indexes at once"
    @assert typeof(color_by) <: Integer
    @assert typeof(design_space) <: Bool
    @assert typeof(pareto) <: Bool
    @assert typeof(max_samples) <: Union{Nothing,Integer}
    @assert typeof(iterations) <: Union{Nothing,Integer}

    if design_space
        arg = :x
        col = :f
        arg_labels = results.opt_ini
        col_labels = results.objectives_functions
    else
        arg = :f
        col = :x
        arg_labels = results.objectives_functions
        col_labels = results.opt_ini
    end

    x = Float64[]
    y = Float64[]
    z = Float64[]
    c = Float64[]

    @series begin
        seriestype --> :scatter
        label --> ""

        # labels
        sol = results.state.convergence[1].population
        if length(indexes) == 1 || length(sol[1].x) == 1
            xlabel --> "Run number"
            ylabel --> pretty_label(arg_labels[indexes[1]])
        end
        if length(indexes) >= 1 && length(sol[1].x) >= 1
            xlabel --> pretty_label(arg_labels[indexes[1]])
        end
        if length(indexes) >= 2 && length(sol[1].x) >= 2
            ylabel --> pretty_label(arg_labels[indexes[2]])
        end
        if length(indexes) == 3 && length(sol[1].x) >= 3
            zlabel --> pretty_label(arg_labels[indexes[3]])
        end

        for (iteration, res) in enumerate(results.state.convergence)
            if iterations !== nothing && iteration ∉ iterations
                continue
            end
            if pareto
                sol = Metaheuristics.get_non_dominated_solutions(res.population)
            else
                sol = res.population
            end

            # data
            if length(indexes) >= 1 && length(sol[1].x) >= 1
                append!(x, (getfield(s, arg)[indexes[1]] for s in sol))
            end
            if length(indexes) >= 2 && length(sol[1].x) >= 2
                append!(y, (getfield(s, arg)[indexes[2]] for s in sol))
            end
            if length(indexes) == 3 && length(sol[1].x) >= 3
                append!(z, (getfield(s, arg)[indexes[3]] for s in sol))
            end
            if color_by == 0
                append!(c, (iteration for s in sol))
            elseif color_by > 0
                append!(c, (getfield(s, col)[color_by] for s in sol))
            end
        end

        # need to go from cost to optimization function
        if !design_space
            if length(indexes) >= 1 && length(sol[1].x) >= 1
                x = results.objectives_functions[indexes[1]].(x)
            end
            if length(indexes) >= 2 && length(sol[1].x) >= 2
                y = results.objectives_functions[indexes[2]].(y)
            end
            if length(indexes) == 3 && length(sol[1].x) >= 3
                z = results.objectives_functions[indexes[3]].(z)
            end
        elseif color_by > 0
            c = results.objectives_functions[color_by].(c)
        end

        # subsample (plotlyjs 3D scatter with large number of points is can be very slow)
        index = collect(1:length(x))
        if max_samples !== nothing
            index = Random.shuffle!(index)[1:max_samples]
            sort!(index)
        end

        # coloring
        if length(c) > 0
            marker_z --> c[index]
            if color_by == 0
                title --> "Colored by iteration"
            else
                title --> "Colored by " * pretty_label(col_labels[color_by])
            end
        end

        # series
        if isempty(y)
            x[index]
        elseif isempty(z)
            x[index], y[index]
        else
            x[index], y[index], z[index]
        end
    end
end

# Everything below is necassary to allow BSON saving/loading of MultiobjectiveOptimizationResults
import IMASDD
import OrderedCollections

function Base.convert(::Type{Vector{<:Parameter}}, x::Vector{Any})
    return Parameter[xx for xx in x]
end

function Base.convert(::Type{Vector{<:ObjectiveFunction}}, x::Vector{Any})
    return ObjectiveFunction[xx for xx in x]
end
