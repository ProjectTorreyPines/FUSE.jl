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
    # inner constructor to register ObjectiveFunction in ObjectiveFunctionsLibrary
    ObjectiveFunction(name::Symbol, units::String, func::Function, target::Float64) = begin
        objf = new(name, units, func, target)
        ObjectiveFunctionsLibrary[objf.name] = objf
        return objf
    end
end

const ObjectiveFunctionsLibrary = Dict{Symbol,ObjectiveFunction}()
function update_ObjectiveFunctionsLibrary!()
    empty!(ObjectiveFunctionsLibrary)
    ObjectiveFunction(:min_levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE, -Inf)
    ObjectiveFunction(:min_log10_levelized_CoE, "log₁₀(\$/kW)", dd -> log10(dd.costing.levelized_CoE), -Inf)
    ObjectiveFunction(:min_capital_cost, "\$B", dd -> dd.costing.cost_direct_capital.cost / 1E3, -Inf)
    ObjectiveFunction(:max_fusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6, Inf)
    ObjectiveFunction(:max_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6, Inf)
    ObjectiveFunction(:max_flattop, "hours", dd -> dd.build.oh.flattop_duration / 3600.0, Inf)
    ObjectiveFunction(:max_log10_flattop, "log₁₀(hours)", dd -> log10(dd.build.oh.flattop_duration / 3600.0), Inf)
    ObjectiveFunction(:min_βn, "", dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal, -Inf)
    return ObjectiveFunctionsLibrary
end
update_ObjectiveFunctionsLibrary!()

"""
    (objf::ObjectiveFunction)(dd::IMAS.dd)

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

function Base.show(io::IO, f::ObjectiveFunction)
    printstyled(io, f.name; bold=true, color=:blue)
    print(io, " →")
    print(io, " $(f.target)")
    print(io, " [$(f.units)]")
end

function Base.show(io::IO, x::MIME"text/plain", objfs::AbstractDict{Symbol,ObjectiveFunction})
    for objf in objfs
        show(io, x, objf)
        println(io, "")
    end
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

const ConstraintFunctionsLibrary = Dict{Symbol,ConstraintFunction}() #s
function update_ConstraintFunctionsLibrary!()
    empty!(ConstraintFunctionsLibrary)
    ConstraintFunction(:target_Beta_n, "", dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal, ==, NaN, 1e-2)
    ConstraintFunction(:target_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6, ==, NaN, 1e-2)
    ConstraintFunction(:steady_state, "log₁₀(hours)", dd -> log10(dd.build.oh.flattop_duration / 3600.0), >, log10(10.0))
    ConstraintFunction(:zero_ohmic, "MA", dd -> abs(sum(integrate(dd.core_profiles.profiles_1d[].grid.area, dd.core_profiles.profiles_1d[].j_ohmic))) / 1E6, ==, 0.0, 1E-2)
    return ConstraintFunctionsLibrary
end
update_ConstraintFunctionsLibrary!()

function (cnst::ConstraintFunction)(dd::IMAS.dd)
    if ===(cnst.operation, ==)
        if cnst.limit === 0.0
            return abs((cnst.func(dd) - cnst.limit))# - cnst.tolerance
        else
            return abs((cnst.func(dd) - cnst.limit) / cnst.limit) - cnst.tolerance
        end
    elseif cnst.operation(1.0, 0.0) # > or >=
        return cnst.limit - cnst.func(dd)
    else # < or <=
        return cnst.func(dd) - cnst.limit
    end
end

function Base.show(io::IO, cnst::ConstraintFunction)
    printstyled(io, cnst.name; bold=true, color=:blue)
    print(io, " $(cnst.operation)")
    print(io, " $(cnst.limit)")
    if ===(cnst.operation, ==)
        if cnst.limit == 0.0
            print(io, " ± $(cnst.tolerance)")
        else
            print(io, " ± $(cnst.tolerance * cnst.limit)")
        end
    end
    print(io, " [$(cnst.units)]")
end

function Base.show(io::IO, x::MIME"text/plain", cnsts::AbstractDict{Symbol,ConstraintFunction})
    for cnst in cnsts
        show(io, x, cnst)
        println(io, "")
    end
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
    constraints_functions::AbstractVector{<:ConstraintFunction},
    save_folder::AbstractString
)
    # update ini based on input optimization vector `x`
    for (optpar, xx) in zip(opt_ini, x)
        if typeof(optpar.value) <: Integer
            optpar.value = Int(round(xx))
        else
            optpar.value = xx
        end
    end
    # run the problem
    try
        if typeof(actor_or_workflow) <: DataType
            actor = actor_or_workflow(init(ini, act), act)
            dd = actor.dd
        else
            dd = actor_or_workflow(ini, act)
        end
        # save simulation data to directory
        if !isempty(save_folder)
            savedir = joinpath(save_folder, "$(Dates.now())__$(getpid())")
            save(savedir, dd, ini, act; freeze=true)
        end
        # evaluate multiple objectives
        return collect(map(f -> nan2inf(f(dd)), objectives_functions)), collect(map(g -> nan2inf(g(dd)), constraints_functions)), Float64[]
    catch e
        # save empty dd and error to directory
        if !isempty(save_folder)
            if typeof(e) <: Exception # somehow sometimes `e` is of type String?
                savedir = joinpath(save_folder, "$(Dates.now())__$(getpid())")
                save(savedir, IMAS.dd(), ini, act, e; freeze=true)
            else
                @warn "typeof(e) in optimization_engine is String: $e"
            end
        end
        # rethrow() # uncomment for debugging purposes
        return Float64[Inf for f in objectives_functions], Float64[Inf for g in constraints_functions], Float64[]
    end
end

"""
    nan2inf(x::Float64)::Float64

Turn NaNs into Inf
"""
function nan2inf(x::Float64)::Float64
    if isnan(x)
        return Inf
    else
        return x
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
    save_folder::AbstractString,
    p::ProgressMeter.Progress)

    # parallel evaluation of a generation
    ProgressMeter.next!(p)
    tmp = Distributed.pmap(x -> optimization_engine(ini, act, actor_or_workflow, x, opt_ini, objectives_functions, constraints_functions, save_folder), [X[k, :] for k in 1:size(X)[1]])
    F = zeros(size(X)[1], length(tmp[1][1]))
    G = zeros(size(X)[1], max(length(tmp[1][2]), 1))
    H = zeros(size(X)[1], max(length(tmp[1][3]), 1))
    for k in 1:size(X)[1]
        f, g, h = tmp[k]
        F[k, :] .= f
        if !isempty(g)
            G[k, :] .= g
        end
        if !isempty(h)
            H[k, :] .= h
        end
    end
    return F, G, H
end