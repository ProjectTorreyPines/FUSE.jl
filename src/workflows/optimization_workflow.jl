function workflow_multiobjective_optimization(func::Function, dd::IMAS.dd, ini::ParametersInit, act::ParametersActor, objectives_functions::AbstractVector{T}=ObjectiveFunction[]; N=10, iterations=N) where {T<:ObjectiveFunction}
    display("Running on $(nprocs()) processes")
    if isempty(objectives_functions)
        error("Must specify objective functions. Available pre-baked functions from ObjectivesFunctionsLibrary:\n  * " * join(keys(ObjectivesFunctionsLibrary), "\n  * "))
    end
    # itentify optimization variables in ini
    opt_ini = opt_parameters(ini)
    for optpar in opt_ini
        display(optpar)
    end
    # optimization boundaries
    bounds = [[optpar.lower for optpar in opt_ini] [optpar.upper for optpar in opt_ini]]'
    # test running function with nominal parameters
    func(dd, ini, act)
    # optimize
    options = Metaheuristics.Options(seed=1, parallel_evaluation=true, store_convergence = true, iterations=iterations)
    algorithm = Metaheuristics.NSGA2(; N, options)
    p = Progress(iterations; desc="Iteration", showspeed=true)
    @time state = Metaheuristics.optimize(X -> optimization_engine(func, dd, ini, act, X, opt_ini, objectives_functions, p), bounds, algorithm)
    return state
end
