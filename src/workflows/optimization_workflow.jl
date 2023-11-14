import Distributed

"""
    workflow_multiobjective_optimization(
        ini::ParametersAllInits,
        act::ParametersAllActors,
        actor_or_workflow::Union{Type{<:AbstractActor},Function},
        objectives_functions::Vector{<:ObjectiveFunction}=ObjectiveFunction[],
        constraints_functions::Vector{<:ConstraintFunction}=ConstraintFunction[];
        exploitation_vs_exploration::Float64=0.0,
        N::Int=10,
        iterations::Int=N,
        continue_state::Union{Missing,Metaheuristics.State}=missing,
        save_folder::AbstractString="optimization_runs",
        save_dd::Bool=true)

Multi-objective optimization of either an `actor(dd, act)` or a `workflow(ini, act)`
"""
function workflow_multiobjective_optimization(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    objectives_functions::Vector{<:ObjectiveFunction}=ObjectiveFunction[],
    constraints_functions::Vector{<:ConstraintFunction}=ConstraintFunction[];
    η_cr::Int=5,
    p_cr::Float64=0.9,
    η_m::Int=5,
    p_m::Float64=1.0,
    algorithm::AbstractString="SPEA2",
    N::Int=10,
    iterations::Int=N,
    continue_state::Union{Missing,Metaheuristics.State}=missing,
    save_folder::AbstractString="optimization_runs",
    save_dd::Bool=true)

    if mod(N, 2) > 0
        error("workflow_multiobjective_optimization population size `N` must be an even number")
    end

    println("Running on $(Distributed.nprocs()-1) worker processes")
    if isempty(objectives_functions)
        error(
            "Must specify objective functions. Available pre-baked functions from ObjectiveFunctionsLibrary:\n  * " *
            join(keys(ObjectiveFunctionsLibrary), "\n  * ")
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
    println("== Constraints ==")
    for cnst in constraints_functions
        println(cnst)
    end

    # scale mutation probability by dimensionality of problem
    D = length(opt_ini)
    p_m = p_m/D

    # optimization boundaries
    bounds = [[float_bounds(optpar)[1] for optpar in opt_ini] [float_bounds(optpar)[2] for optpar in opt_ini]]'

    # # test running function once with nominal parameters useful to catch bugs quickly.
    # # Use Distributed.@everywhere to trigger compilation on all worker nodes.
    # if typeof(actor_or_workflow) <: Function
    #     actor_or_workflow(ini, act)
    # else
    #     actor_or_workflow(init(ini, act), act)
    # end

    save_folder = abspath(save_folder)
    if !isempty(save_folder)
        mkpath(save_folder)
    end

    # optimize
    options = Metaheuristics.Options(; iterations, parallel_evaluation=true, store_convergence=true, seed=1, f_calls_limit=1E9, g_calls_limit=1E9, h_calls_limit=1E9)

    println()
    println("== Algorithm settings ==")
    println("Algorithm name: "*algorithm)
    println("Crossover distribution index: η_cr = "*string(η_cr))
    println("Crossover probability: p_cr = "*string(p_cr))
    println("Mutation distribution index: η_m = "*string(η_m))
    println("Mutation probability: p_m = "*string(p_m))
    println()


    if algorithm == "SPEA2"
        algorithm_obj = Metaheuristics.SPEA2(; N, η_cr, p_cr, η_m, p_m, options)
    elseif algorithm == "NSGA2"
        algorithm_obj = Metaheuristics.NSGA2(; N, η_cr, p_cr, η_m, p_m, options) 
    elseif algorithm == "CCMO"
        algorithm_obj = Metaheuristics.CCMO(Metaheuristics.NSGA2(; N, η_cr, p_cr, η_m, p_m, options)) 
    end
    if continue_state !== missing
        println("Restarting simulation")
        algorithm_obj.status = continue_state
    end
    flush(stdout)

    p = ProgressMeter.Progress(iterations; desc="Iteration", showspeed=true)
    @time state =
        Metaheuristics.optimize(X -> optimization_engine(ini, act, actor_or_workflow, X, objectives_functions, constraints_functions, save_folder, save_dd, p), bounds, algorithm)
    display(state)

    return state
end

function is_dominated(sol_a::Vector{T}, sol_b::Vector{T}) where {T}
    return all(sol_b .<= sol_a) .&& any(sol_b .< sol_a)
end

"""
    pareto_front(solutions::Vector{Vector{T}}) where T

returns indexes of solutions that form the pareto front
"""
function pareto_front(solutions::Vector{Vector{T}}) where {T}
    pareto = Int[]
    for i in eachindex(solutions)
        is_dominated_by_any = false
        for j in eachindex(solutions)
            if i != j && is_dominated(solutions[i], solutions[j])
                is_dominated_by_any = true
                break
            end
        end
        if !is_dominated_by_any
            push!(pareto, i)
        end
    end
    return pareto
end
