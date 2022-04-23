import Metaheuristics

function optimization_engine(func::Function, dd::IMAS.dd, ini::InitParameters, act::ActorParameters, x::AbstractVector, objectives_functions::AbstractVector{T}) where {T<:Function}
    # update ini based on input optimization vector `x`
    for (optpar, xx) in zip(opt_parameters(ini), x)
        if typeof(optpar.value) <: Integer
            optpar.value = Int(round(xx))
        else
            optpar.value = xx
        end
    end
    # call the problem
    func(dd, ini, act)
    # evaluate multiple objectives
    return collect(map(f -> f(dd), objectives_functions)), x * 0, x * 0
end

function optimization_engine(func::Function, dd::IMAS.dd, ini::InitParameters, act::ActorParameters, X::AbstractMatrix, objectives_functions::AbstractVector{T}) where {T<:Function}
    display("Running on $(nprocs()) processes")
    # parallel evaluation of a generation
    tmp = pmap(x -> optimization_engine(func, dd, ini, act, x, objectives_functions), [@view X[k, :] for k in 1:size(X)[1]])
    F = zeros(size(X)[1], length(tmp[1][1]))
    G = similar(X)
    H = similar(X)
    for k in 1:size(X)[1]
        F[k, :], G[k, :], H[k, :] = tmp[k][1], tmp[k][2], tmp[k][3]
    end
    return F, G, H
end

function optimization_workflow(func::Function, dd::IMAS.dd, ini::InitParameters, act::ActorParameters, objectives_functions::AbstractVector{T}=Function[]; N=10) where {T<:Function}
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
    # optimize
    options = Metaheuristics.Options(parallel_evaluation=true, f_calls_limit=10 * N)
    algorithm = Metaheuristics.NSGA2(N=N, options=options)
    @time results = Metaheuristics.optimize(X -> optimization_engine(func, dd, ini, act, X, objectives_functions), bounds, algorithm)
    return results
end

ObjectivesFunctionsLibrary = Dict{Symbol,Function}()
ObjectivesFunctionsLibrary[:min_cost] = dd -> dd.costing.cost
ObjectivesFunctionsLibrary[:max_fusion] = dd -> -IMAS.fusion_power(dd.core_profiles.profiles_1d[])
ObjectivesFunctionsLibrary[:min_ohmic_current] = dd -> @ddtime dd.summary.global_quantities.current_ohm.value