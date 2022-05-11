import Distributed
import OrderedCollections # necassary to allow BSON saving of MultiobjectiveOptimizationResults

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
    continue_results::Union{Missing,MultiobjectiveOptimizationResults}=missing)

    if mod(N, 2) > 0
        error("workflow_multiobjective_optimization population size `N` must be an even number")
    end

    println("Running on $(nprocs()) processes")
    if isempty(objectives_functions)
        error("Must specify objective functions. Available pre-baked functions from ObjectivesFunctionsLibrary:\n  * " * join(keys(ObjectivesFunctionsLibrary), "\n  * "))
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
        actor_or_workflow(FUSE.init(ini, act), act)
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

function pretty_label(objective_function::ObjectiveFunction, units="")
    txt = split(string(objective_function.name), "_")[2]
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

@recipe function plot_MultiobjectiveOptimizationResults(results::MultiobjectiveOptimizationResults, indexes=[1, 2, 3]; design_space::Bool=false, pareto=true)

    if design_space
        arg = :x
        labels = results.opt_ini
    else
        arg = :f
        labels = results.objectives_functions
    end

    x = []
    y = []
    z = []
    c = []

    @series begin
        seriestype --> :scatter
        label --> ""
        for (generation, res) in enumerate(results.state.convergence)
            if pareto
                sol = Metaheuristics.get_non_dominated_solutions(res.population)
            else
                sol = res.population
            end

            if length(indexes) == 1 || length(sol[1].x) == 1
                xlabel --> "Run number"
                ylabel --> pretty_label(labels[indexes[1]])
                append!(x, (getfield(s,arg)[indexes[1]] for s in sol))
            elseif length(indexes) == 2 || length(sol[1].x) == 2
                append!(x, (getfield(s,arg)[indexes[1]] for s in sol))
                append!(y, (getfield(s,arg)[indexes[2]] for s in sol))
                xlabel --> pretty_label(labels[indexes[1]])
                ylabel --> pretty_label(labels[indexes[2]])
            elseif length(indexes) == 3 || length(sol[1].x) == 3
                append!(x, (getfield(s,arg)[indexes[1]] for s in sol))
                append!(y, (getfield(s,arg)[indexes[2]] for s in sol))
                append!(z, (getfield(s,arg)[indexes[3]] for s in sol))
                xlabel --> pretty_label(labels[indexes[1]])
                ylabel --> pretty_label(labels[indexes[2]])
                zlabel --> pretty_label(labels[indexes[3]])
            end

            #append!(c, (s.x[3] for s in sol))
            append!(c, (generation for s in sol))
        end

        if length(c) > 0
            marker_z --> c
        end

        if isempty(y)
            x
        elseif isempty(z)
            x, y
        else
            x, y, z
        end
    end

end