import Metaheuristics

function optimization_engine(func::Function, dd::IMAS.dd, ini::InitParameters, act::ActorParameters, opt_vector, x::AbstractVector)
    for (optpar, xx) in zip(opt_vector, x)
        if typeof(optpar.value) <: Integer
            optpar.value = Int(round(xx))
        else
            optpar.value = xx
        end
    end
    func(dd, ini, act)
    return [abs(dd.core_profiles.profiles_1d[].electrons.temperature[1]-30e3)/dd.core_profiles.profiles_1d[].electrons.temperature[1],
            abs(dd.core_profiles.profiles_1d[].electrons.temperature[1]-29e3)/dd.core_profiles.profiles_1d[].electrons.temperature[1],], x * 0, x * 0
end

function optimization_engine(func::Function, dd::IMAS.dd, ini::InitParameters, act::ActorParameters, opt_vector, X::AbstractMatrix)
    display("Running on $(nprocs()) processes")
    tmp = pmap(x -> optimization_engine(func, dd, ini, act, opt_vector, x), [@view X[k, :] for k in 1:size(X)[1]])
    F = zeros(size(X)[1], length(tmp[1][1]))
    G = similar(X)
    H = similar(X)
    for k in 1:size(X)[1]
        F[k, :], G[k, :], H[k, :] = tmp[k][1], tmp[k][2], tmp[k][3]
    end
    return F, G, H
end

function optimization_workflow(func, dd, ini, act)
    opt_vector = opt_parameters(ini)
    for optpar in opt_vector
        display(optpar)
    end
    bounds = [[optpar.lower for optpar in opt_vector] [optpar.upper for optpar in opt_vector]]'
    options = Metaheuristics.Options(parallel_evaluation=true,f_calls_limit = 10)
    algorithm = Metaheuristics.NSGA2(N=10, options=options)
    @time results = Metaheuristics.optimize(X -> optimization_engine(func, dd, ini, act, opt_vector, X), bounds, algorithm)
    return results
end
