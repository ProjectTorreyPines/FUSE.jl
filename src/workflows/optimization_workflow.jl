import Distributed

mutable struct MultiobjectiveOptimizationResults
    workflow::Union{DataType, Function}
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
    actor_or_workflow::Union{DataType, Function},
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
