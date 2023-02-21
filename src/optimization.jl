import Metaheuristics
import ProgressMeter
import Distributed
import Dates
ProgressMeter.ijulia_behavior(:clear)

# ==================== #
# objectives functions #
# ==================== #
mutable struct ObjectiveFunction
    name::Symbol
    units::String
    func::Function
    target::Float64
    # inner constructor to register ObjectiveFunction in ObjectivesFunctionsLibrary
    ObjectiveFunction(name::Symbol, units::String, func::Function, target::Float64) = begin
        objf = new(name, units, func, target)
        ObjectivesFunctionsLibrary[objf.name] = objf
        return objf
    end
end

"""
    (objf::ObjectiveFunction)(x::Float64)

From real domain to objective domain (Metaheuristics will always minimize)
"""
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

"""
    (objf::ObjectiveFunction)(x::Float64)

From objective domain to real domain
"""
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
ObjectiveFunction(:min_levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE, -Inf)
ObjectiveFunction(:min_log10_levelized_CoE, "log₁₀(\$/kW)", dd -> log10(dd.costing.levelized_CoE), -Inf)
ObjectiveFunction(:max_fusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6, Inf)
ObjectiveFunction(:max_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6, Inf)
ObjectiveFunction(:max_flattop, "hours", dd -> dd.build.oh.flattop_duration / 3600.0, Inf)
ObjectiveFunction(:max_log10_flattop, "log₁₀(hours)", dd -> log10(dd.build.oh.flattop_duration / 3600.0), Inf)

function Base.show(io::IO, f::ObjectiveFunction)
    printstyled(io, f.name; bold=true, color=:blue)
    print(io, " →")
    print(io, " $(f.target)")
    print(io, " [$(f.units)]")
end

# ==================== #
# constraint functions #
# ==================== #
mutable struct ConstraintFunction
    name::Symbol
    units::String
    func::Function
    operation::Function
    limit::Float64
    tolerance::Float64
    # inner constructor to register ConstraintFunction in ConstraintsFunctionsLibrary
    ConstraintFunction(name::Symbol, units::String, func::Function, operation::Function, limit::Float64, tolerance::Float64) = begin
        @assert ===(operation, ==) "tolerance specification only used for == constraint"
        cnst = new(name, units, func, operation, limit, tolerance)
        ConstraintFunctionsLibrary[cnst.name] = cnst
        return cnst
    end
    ConstraintFunction(name::Symbol, units::String, func::Function, operation::Function, limit::Float64) = begin
        @assert !==(operation, ==) "Must specify tolerance of == constraint"
        cnst = new(name, units, func, operation, limit, 0.0)
        ConstraintFunctionsLibrary[cnst.name] = cnst
        return cnst
    end
end

function (cnst::ConstraintFunction)(dd::IMAS.dd)
    if ===(cnst.operation, ==)
        return (cnst.func(dd) - cnst.limit)^2 - (cnst.limit * cnst.tolerance)^2
    elseif cnst.operation(1.0, 0.0) # > or >=
        return  - (cnst.func(dd) - cnst.limit)
    else # < or <=
        return  cnst.func(dd) - cnst.limit
    end
end

const ConstraintFunctionsLibrary = Dict{Symbol,ConstraintFunction}() #s
ConstraintFunction(:target_Beta_n, "", dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal, ==, NaN, 1e-2)
ConstraintFunction(:target_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6, ==, NaN, 1e-2)
ConstraintFunction(:steady_state, "hours", dd -> dd.build.oh.flattop_duration / 3600.0, >, 10.0)

function Base.show(io::IO, f::ConstraintFunction)
    printstyled(io, f.name; bold=true, color=:blue)
    print(io, " $(f.operation)")
    print(io, " $(f.limit)")
    if ===(f.operation, ==)
        print(io, " ± $(f.tolerance * f.limit)")
    end
    print(io, " [$(f.units)]")
end

# =================== #
# Optimization engine #
# =================== #
function optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{DataType,Function},
    x::AbstractVector,
    opt_ini::Vector{<:AbstractParameter},
    objectives_functions::AbstractVector{<:ObjectiveFunction},
    constraints_functions::AbstractVector{<:ConstraintFunction}
    )
    # update ini based on input optimization vector `x`
    for (optpar, xx) in zip(opt_ini, x)
        if typeof(optpar.value) <: Integer
            optpar.value = Int(round(xx))
        else
            optpar.value = xx
        end
    end
    topdirname = "optimization_runs"
    mkpath(topdirname)
    # run the problem
    try
        if typeof(actor_or_workflow) <: DataType
            actor = actor_or_workflow(init(ini, act), act)
            dd = actor.dd
        else
            dd = actor_or_workflow(ini, act)
        end
        # save simulation data to directory
        savedir = joinpath(topdirname, "$(Dates.now())__$(getpid())")
        save(dd, ini, act, savedir; freeze=true)
        # evaluate multiple objectives
        return collect(map(f -> f(dd), objectives_functions)), Float64[], collect(map(h -> h(dd), constraints_functions))
    catch e
        # save empty dd and error to directory
        savedir = joinpath(topdirname, "$(Dates.now())__$(getpid())")
        save(IMAS.dd(), ini, act, savedir; freeze=true)
        open(joinpath(savedir, "error.txt"), "w") do file
            showerror(file, e, catch_backtrace())
        end
        # rethrow() # uncomment for debugging purposes
        return Float64[Inf for f in objectives_functions], Float64[], Float64[Inf for h in constraints_functions]
    end
end

function optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{DataType,Function},
    X::AbstractMatrix,
    opt_ini::Vector{<:AbstractParameter},
    objectives_functions::AbstractVector{<:ObjectiveFunction},
    constraints_functions::AbstractVector{<:ConstraintFunction},
    p::ProgressMeter.Progress)

    # parallel evaluation of a generation
    ProgressMeter.next!(p)
    tmp = Distributed.pmap(x -> optimization_engine(ini, act, actor_or_workflow, x, opt_ini, objectives_functions, constraints_functions), [X[k, :] for k in 1:size(X)[1]])
    F = zeros(size(X)[1], length(objectives_functions))
    G = zeros(size(X)[1], 0)
    H = zeros(size(X)[1], length(constraints_functions))
    for k in 1:size(X)[1]
        F[k, :] .= tmp[k][1]
        G[k, :] .= tmp[k][2]
        H[k, :] .= tmp[k][3]
    end
    return F, G, H
end
