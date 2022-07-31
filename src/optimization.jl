import Metaheuristics
using ProgressMeter
import Distributed
ProgressMeter.ijulia_behavior(:clear)

mutable struct ObjectiveFunction
    name::Symbol
    units::String
    func::Function
    target::Float64
    # inner constructor to register ObjectiveFunction in ObjectivesFunctionsLibrary
    ObjectiveFunction(name, units, func, target) = begin
        objf = new(name, units, func, target)
        ObjectivesFunctionsLibrary[objf.name] = objf
        return objf
    end
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

const ObjectivesFunctionsLibrary = Dict{Symbol,ObjectiveFunction}()
ObjectiveFunction(:min_cost, "\$M", dd -> dd.costing.cost, -Inf)
ObjectiveFunction(:max_fusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6, Inf)
ObjectiveFunction(:max_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6, Inf)
ObjectiveFunction(:max_flattop, "hours", dd -> dd.build.oh.flattop_duration / 3600, Inf)
ObjectiveFunction(:max_log10_flattop, "log₁₀(hours)", dd -> log10(dd.build.oh.flattop_duration / 3600), Inf)

function Base.show(io::IO, f::ObjectiveFunction)
    printstyled(io, f.name; bold=true, color=:blue)
    print(io, " [$(f.units)]")
end

function optimization_engine(ini::ParametersAllInits, act::ParametersAllActors, actor_or_workflow::Union{DataType,Function}, x::AbstractVector, opt_ini, objectives_functions::AbstractVector{<:ObjectiveFunction})
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
        if typeof(actor_or_workflow) <: DataType
            actor = actor_or_workflow(init(ini, act), act)
            dd = actor.dd
        else
            dd = actor_or_workflow(ini, act)
        end
        # evaluate multiple objectives
        return collect(map(f -> f(dd), objectives_functions)), x * 0, x * 0
    catch
        # rethrow() 
        return [Inf for f in objectives_functions], x * 0, x * 0
    end
end

function optimization_engine(ini::ParametersAllInits, act::ParametersAllActors, actor_or_workflow::Union{DataType,Function}, X::AbstractMatrix, opt_ini, objectives_functions::AbstractVector{<:ObjectiveFunction}, p)
    # parallel evaluation of a generation
    ProgressMeter.next!(p)
    tmp = Distributed.pmap(x -> optimization_engine(ini, act, actor_or_workflow, x, opt_ini, objectives_functions), [X[k, :] for k in 1:size(X)[1]])
    F = zeros(size(X)[1], length(objectives_functions))
    G = similar(X)
    H = similar(X)
    for k in 1:size(X)[1]
        F[k, :] .= tmp[k][1]
        G[k, :] .= tmp[k][2]
        H[k, :] .= tmp[k][3]
    end
    return F, G, H
end
