import Metaheuristics
using ProgressMeter
ProgressMeter.ijulia_behavior(:clear)

mutable struct ObjectiveFunction
    name::Symbol
    func::Function
    target::Float64
    registered::Bool
end

function ObjectiveFunction(name::Symbol, func::Function, target::Float64)
    objf = ObjectiveFunction(name, func, target, true)
    ObjectivesFunctionsLibrary[objf.name] = objf
end

function (objf::ObjectiveFunction)(dd::IMAS.dd)
    if isinf(objf.target)
        if objf.target < 0
            return objf.func(dd)
        else
            return -objf.func(dd)
        end
    elseif objf.target == 0.0
        return abs(objf.func(dd))
    else
        return abs(objf.func(dd) - objf.target) / objf.target
    end
end

function (objf::ObjectiveFunction)(x::Float64)
    if isinf(objf.target)
        if objf.target < 0
            return x
        else
            return -x
        end
    elseif objf.target == 0.0
        return x
    else
        return x * objf.target + objf.target
    end
end

function optimization_engine(func::Function, dd::IMAS.dd, ini::InitParameters, act::ActorParameters, x::AbstractVector, opt_ini, objectives_functions::AbstractVector{T}) where {T<:ObjectiveFunction}
    # update ini based on input optimization vector `x`
    for (optpar, xx) in zip(opt_ini, x)
        if typeof(optpar.value) <: Integer
            optpar.value = Int(round(xx))
        else
            optpar.value = xx
        end
    end
    # call the problem
    try
        func(dd, ini, act)
    catch
        #return [Inf for f in objectives_functions], x * 0, x * 0
        rethrow()
    end
    # evaluate multiple objectives
    return collect(map(f -> f(dd), objectives_functions)), x * 0, x * 0
end

function optimization_engine(func::Function, dd::IMAS.dd, ini::InitParameters, act::ActorParameters, X::AbstractMatrix, opt_ini, objectives_functions::AbstractVector{T}, p) where {T<:ObjectiveFunction}
    # parallel evaluation of a generation
    ProgressMeter.next!(p)
    tmp = pmap(x -> optimization_engine(func, dd, ini, act, x, opt_ini, objectives_functions), [X[k, :] for k in 1:size(X)[1]])
    F = zeros(size(X)[1], length(tmp[1][1]))
    G = similar(X)
    H = similar(X)
    for k in 1:size(X)[1]
        F[k, :], G[k, :], H[k, :] = tmp[k][1], tmp[k][2], tmp[k][3]
    end
    return F, G, H
end

function optimization_workflow(func::Function, dd::IMAS.dd, ini::InitParameters, act::ActorParameters, objectives_functions::AbstractVector{T}=ObjectiveFunction[]; N=10, iterations=N) where {T<:ObjectiveFunction}
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
    options = Metaheuristics.Options(parallel_evaluation=true, iterations=iterations)
    algorithm = Metaheuristics.NSGA2(N=N, options=options)
    p = Progress(iterations; desc="Iteration", showspeed=true)
    @time state = Metaheuristics.optimize(X -> optimization_engine(func, dd, ini, act, X, opt_ini, objectives_functions, p), bounds, algorithm)
    return state
end

ObjectivesFunctionsLibrary = Dict{Symbol,ObjectiveFunction}()
ObjectiveFunction(:min_cost, dd -> dd.costing.cost, -Inf)
ObjectiveFunction(:max_fusion, dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]), Inf)
ObjectiveFunction(:max_flattop_duration, dd -> bd.oh.flattop_duration, Inf)
