import Distributed
import Serialization

"""
    workflow_multiobjective_optimization(
        ini::ParametersAllInits,
        act::ParametersAllActors,
        actor_or_workflow::Union{Type{<:AbstractActor},Function},
        objective_functions::Vector{<:IMAS.ObjectiveFunction}=ObjectiveFunction[],
        constraint_functions::Vector{<:ConstraintFunction}=ConstraintFunction[];
        exploitation_vs_exploration::Float64=0.0,
        N::Int=10,
        iterations::Int=N,
        continue_state::Union{Nothing,Metaheuristics.State}=missing,
        save_folder::AbstractString="optimization_runs",
        save_dd::Bool=true)

Multi-objective optimization of either an `actor(dd, act)` or a `workflow(ini, act)`
"""
function workflow_multiobjective_optimization(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    objective_functions::Vector{<:IMAS.ObjectiveFunction}=IMAS.ObjectiveFunction[],
    constraint_functions::Vector{<:IMAS.ConstraintFunction}=IMAS.ConstraintFunction[];
    exploitation_vs_exploration::Float64=0.0,
    N::Int=10,
    iterations::Int=N,
    continue_state::Union{Nothing,Metaheuristics.State}=nothing,
    save_folder::AbstractString="optimization_runs",
    save_dd::Bool=true,
    generation_offset::Int=0)

    if mod(N, 2) > 0
        error("workflow_multiobjective_optimization population size `N` must be an even number")
    end

    println("Running on $(Distributed.nprocs()-1) worker processes")
    if isempty(objective_functions)
        error(
            "Must specify objective functions. Available pre-baked functions from ObjectiveFunctionsLibrary:\n  * " *
            join(keys(IMAS.ObjectiveFunctionsLibrary), "\n  * ")
        )
    end

    # itentify optimization variables in ini
    opt_ini = opt_parameters(ini)
    if generation_offset == 0
        println("== Actuators ==")
        for optpar in opt_ini
            println(optpar)
        end
        println()
        println("== Objectives ==")
        for objf in objective_functions
            println(objf)
        end
        println()
        println("== Constraints ==")
        for cnst in constraint_functions
            println(cnst)
        end
    end

    # optimization floating point boundaries
    bounds = float_bounds(opt_ini)

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

    if length(objective_functions) == 1
        D = length(opt_ini)
        K = Int(floor(N / D))
        algorithm = Metaheuristics.ECA(; N, K, options)

    else
        # set algorithm parameters depending on exploitation_vs_exploration index
        # crossover distribution index
        η_cr = Int(round(IMAS.interp1d([0.0, 1.0, 2.0], [20.0, 30.0, 40.0], :cubic).(exploitation_vs_exploration)))
        # crossover probability
        p_cr = IMAS.interp1d([0.0, 1.0, 2.0], [0.9, 0.6, 0.5], :cubic).(exploitation_vs_exploration)
        # mutation distribution index
        η_m = Int(round(IMAS.interp1d([0.0, 1.0, 2.0], [20.0, 30.0, 50.0], :cubic).(exploitation_vs_exploration)))
        # mutation probability
        p_m = IMAS.interp1d([0.0, 1.0, 2.0], [1.0, 2.0, 4.0], :cubic).(exploitation_vs_exploration)

        # algorithm = Metaheuristics.NSGA2(; N, options) # converges to one point and does not cover well the pareto front
        # algorithm = Metaheuristics.SMS_EMOA(; N, options) # does not converge
        # algorithm = Metaheuristics.CCMO(Metaheuristics.NSGA2(; N, options); options) # not better than SPEA2
        algorithm = Metaheuristics.SPEA2(; N, η_cr, p_cr, η_m, p_m, options) # converges and covers well the pareto front!
    end

    if continue_state !== nothing
        println("Restarting simulation")
        algorithm.status = continue_state
    end

    flush(stdout)
    ProgressMeter.ijulia_behavior(:clear)
    p = ProgressMeter.Progress(iterations; desc="Iteration", showspeed=true)
    @time state =
        Metaheuristics.optimize(X -> optimization_engine(ini, act, actor_or_workflow, X, objective_functions, constraint_functions, save_folder, save_dd, p, generation_offset), bounds, algorithm)
    if generation_offset == 0
        display(state)
    end

    if !isempty(save_folder)
        save_optimization(joinpath(save_folder, "results.jls"), state, ini, act, objective_functions, constraint_functions)
    end

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

"""
    save_optimization(
        filename::AbstractString,
        state::Metaheuristics.State,
        ini::ParametersAllInits,
        act::ParametersAllActors,
        objectives_functions::Vector{<:IMAS.ObjectiveFunction},
        constraints_functions::Vector{<:IMAS.ConstraintFunction})

Save Metaheuristics.State to file
"""
function save_optimization(
    filename::AbstractString,
    state::Metaheuristics.State,
    ini::ParametersAllInits,
    act::ParametersAllActors,
    objectives_functions::Vector{<:IMAS.ObjectiveFunction},
    constraints_functions::Vector{<:IMAS.ConstraintFunction}
)
    data = Dict("state" => state, "ini" => ini, "act" => act, "objective_functions" => objectives_functions, "constraint_functions" => constraints_functions)
    open(filename, "w") do io
        return Serialization.serialize(io, data)
    end
end

"""
    load_optimization(filename::AbstractString)

Load Metaheuristics.State from file, returns named tuple
"""
function load_optimization(filename::AbstractString)
    data = open(filename, "r") do io
        return Serialization.deserialize(io)
    end
    return (state=data["state"], ini=data["ini"], act=data["act"], objectives_functions=data["objective_functions"], constraints_functions=data["constraint_functions"])
end
