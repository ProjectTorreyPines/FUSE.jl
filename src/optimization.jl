using Metaheuristics: Metaheuristics
using ProgressMeter: ProgressMeter
using Distributed: Distributed
using Dates: Dates
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
    #! format: off
    ObjectiveFunction(:min_levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE, -Inf)
    ObjectiveFunction(:min_log10_levelized_CoE, "log₁₀(\$/kW)", dd -> log10(dd.costing.levelized_CoE), -Inf)
    ObjectiveFunction(:min_capital_cost, "\$B", dd -> dd.costing.cost_direct_capital.cost / 1E3, -Inf)
    ObjectiveFunction(:max_fusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6, Inf)
    ObjectiveFunction(:max_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6, Inf)
    ObjectiveFunction(:req_power_electric_net, "ΔMW", dd -> abs(@ddtime(dd.balance_of_plant.power_electric_net) - dd.requirements.power_electric_net) / 1E6, 0.0)
    ObjectiveFunction(:max_flattop, "hours", dd -> dd.build.oh.flattop_duration / 3600.0, Inf)
    ObjectiveFunction(:req_flattop, "Δhours", dd -> abs(dd.build.oh.flattop_duration - dd.requirements.flattop_duration) / 3600.0, 0.0)
    ObjectiveFunction(:max_log10_flattop, "log₁₀(hours)", dd -> log10(dd.build.oh.flattop_duration / 3600.0), Inf)
    ObjectiveFunction(:min_βn, "", dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal, -Inf)
    ObjectiveFunction(:min_R0, "m", dd -> dd.equilibrium.time_slice[].boundary.geometric_axis.r, -Inf)
    #! format: on
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
    return print(io, " [$(f.units)]")
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
    #! format: off
    ConstraintFunction(:required_power_electric_net, "%", dd -> abs(@ddtime(dd.balance_of_plant.power_electric_net) - dd.requirements.power_electric_net) / dd.requirements.power_electric_net, ==, 0.0, 0.01) # relative tolerance
    ConstraintFunction(:min_required_power_electric_net, "%", dd -> (@ddtime(dd.balance_of_plant.power_electric_net) - dd.requirements.power_electric_net) / dd.requirements.power_electric_net, >, 0.0)
    ConstraintFunction(:required_flattop, "%", dd -> abs(dd.build.oh.flattop_duration - dd.requirements.flattop_duration) / dd.requirements.flattop_duration, ==, 0.0, 0.01) # relative tolerance
    ConstraintFunction(:min_required_flattop, "%", dd -> (dd.build.oh.flattop_duration - dd.requirements.flattop_duration) / dd.requirements.flattop_duration, >, 0.0)
    ConstraintFunction(:min_required_B0, "%", dd -> (abs(prod(IMAS.build_max_R0_B0(dd.build))) - dd.equilibrium.vacuum_toroidal_field.r0 * maximum(abs, dd.equilibrium.vacuum_toroidal_field.b0)) / (dd.equilibrium.vacuum_toroidal_field.r0 * maximum(abs, dd.equilibrium.vacuum_toroidal_field.b0)), >, 0.0)
    ConstraintFunction(:zero_ohmic, "MA", dd -> abs(sum(integrate(dd.core_profiles.profiles_1d[].grid.area, dd.core_profiles.profiles_1d[].j_ohmic))) / 1E6, ==, 0.0, 0.1) # absolute tolerance
    ConstraintFunction(:max_ne_peaking, "%", dd -> ((@ddtime(dd.summary.local.magnetic_axis.n_e.value) / @ddtime(dd.summary.volume_average.n_e.value)) - dd.requirements.ne_peaking) / dd.requirements.ne_peaking, <, 0.0)
    ConstraintFunction(:min_lh_power_threshold, "%", dd -> (IMAS.power_sol(dd) / dd.requirements.lh_power_threshold_fraction - IMAS.scaling_L_to_H_power(dd)) / IMAS.scaling_L_to_H_power(dd), >, 0.0)
    ConstraintFunction(:max_ωpe_ωce, "%", dd -> IMAS.ω_pe(@ddtime(dd.summary.local.magnetic_axis.n_e.value)) / IMAS.ω_ce(@ddtime(dd.equilibrium.vacuum_toroidal_field.b0)), <, 1.0)
    ConstraintFunction(:max_qpol_omp, "%", dd -> (IMAS.q_pol_omp_eich(dd) - dd.requirements.q_pol_omp) / dd.requirements.q_pol_omp, <, 0.0)
    ConstraintFunction(:max_tf_j, "%", dd -> dd.build.tf.critical_j / dd.build.tf.max_j - dd.requirements.coil_j_margin, >, 0.0)
    ConstraintFunction(:max_oh_j, "%", dd -> dd.build.oh.critical_j / dd.build.oh.max_j - dd.requirements.coil_j_margin, >, 0.0)
    ConstraintFunction(:max_pl_stress, "%", dd -> ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl) ? 0.0 : dd.solid_mechanics.center_stack.properties.yield_strength.pl / maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl) - dd.requirements.coil_stress_margin, >, 0.0)
    ConstraintFunction(:max_tf_stress, "%", dd -> dd.solid_mechanics.center_stack.properties.yield_strength.tf / maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf) - dd.requirements.coil_stress_margin, >, 0.0)
    ConstraintFunction(:max_oh_stress, "%", dd -> dd.solid_mechanics.center_stack.properties.yield_strength.oh / maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh) - dd.requirements.coil_stress_margin, >, 0.0)
    ConstraintFunction(:max_hds03, "%", dd -> (@ddtime(dd.summary.global_quantities.tau_energy.value)/IMAS.tau_e_ds03(dd) - dd.requirements.hds03)/dd.requirements.hds03, <, 0.0)
    ConstraintFunction(:min_q95, "%", dd -> (dd.equilibrium.time_slice[].global_quantities.q_95 - dd.requirements.q95)/dd.requirements.q95, >, 0.0)
    #! format: on
    return ConstraintFunctionsLibrary
end
update_ConstraintFunctionsLibrary!()

function (cnst::ConstraintFunction)(dd::IMAS.dd)
    if ===(cnst.operation, ==)
        if cnst.limit === 0.0
            return abs((cnst.func(dd) - cnst.limit)) - cnst.tolerance
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
    return print(io, " [$(cnst.units)]")
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
"""
    optimization_engine(
        ini::ParametersAllInits,
        act::ParametersAllActors,
        actor_or_workflow::Union{Type{<:AbstractActor},Function},
        x::AbstractVector,
        objective_functions::AbstractVector{<:ObjectiveFunction},
        constraint_functions::AbstractVector{<:ConstraintFunction},
        save_folder::AbstractString,
        generation::Int,
        save_dd::Bool=true)

NOTE: This function is run by the worker nodes
"""
function optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    x::AbstractVector,
    objective_functions::AbstractVector{<:ObjectiveFunction},
    constraint_functions::AbstractVector{<:ConstraintFunction},
    save_folder::AbstractString,
    generation::Int,
    save_dd::Bool=true)

    # create working directory
    original_dir = pwd()
    savedir = abspath(joinpath(save_folder, "$(generation)__$(Dates.now())__$(getpid())"))
    mkdir(savedir)
    cd(savedir)

    # Redirect stdout and stderr to the file
    original_stdout = stdout  # Save the original stdout
    original_stderr = stderr  # Save the original stderr
    file_log = open("log.txt", "w")

    try
        redirect_stdout(file_log)
        redirect_stderr(file_log)

        # deepcopy ini/act to avoid changes
        ini = deepcopy(ini)
        act = deepcopy(act)

        # update ini based on input optimization vector `x`
        @show x
        parameters_from_opt!(ini, x)
        @show x

        # attempt to release memory
        malloc_trim_if_glibc()

        # run the problem
        if typeof(actor_or_workflow) <: Function
            dd = actor_or_workflow(ini, act)
        else
            dd = init(ini, act)
            actor = actor_or_workflow(dd, act)
            dd = actor.dd
        end

        # save simulation data to directory
        save(savedir, save_dd ? dd : nothing, ini, act; timer=true, freeze=false, overwrite_files=true)

        # evaluate multiple objectives
        ff = collect(map(f -> nan2inf(f(dd)), objective_functions))
        gg = collect(map(g -> nan2inf(g(dd)), constraint_functions))
        hh = Float64[]

        println("finished")
        @show ff
        @show gg
        @show hh

        return ff, gg, hh

    catch e
        # save empty dd and error to directory
        save(savedir, nothing, ini, act, e; timer=true, freeze=false, overwrite_files=true)
        
        # rethrow(e) # uncomment for debugging purposes
        
        ff = Float64[Inf for f in objective_functions]
        gg = Float64[Inf for g in constraint_functions]
        hh = Float64[]

        println("failed")
        @show ff
        @show gg
        @show hh

        return ff, gg, hh

    finally
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        cd(original_dir)
        close(file_log)
    end
end

function _optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    x::AbstractVector,
    objective_functions::AbstractVector{<:ObjectiveFunction},
    constraint_functions::AbstractVector{<:ConstraintFunction},
    save_folder::AbstractString,
    generation::Int,
    save_dd::Bool=true)

    tmp = optimization_engine(ini, act, actor_or_workflow, x, objective_functions, constraint_functions, save_folder, generation, save_dd)

    GC.gc()
    return tmp
end

"""
    optimization_engine(
        ini::ParametersAllInits,
        act::ParametersAllActors,
        actor_or_workflow::Union{Type{<:AbstractActor},Function},
        X::AbstractMatrix,
        objective_functions::AbstractVector{<:ObjectiveFunction},
        constraint_functions::AbstractVector{<:ConstraintFunction},
        save_folder::AbstractString,
        save_dd::Bool,
        p::ProgressMeter.Progress)

NOTE: this function is run by the master process
"""
function optimization_engine(
    ini::ParametersAllInits,
    act::ParametersAllActors,
    actor_or_workflow::Union{Type{<:AbstractActor},Function},
    X::AbstractMatrix,
    objective_functions::AbstractVector{<:ObjectiveFunction},
    constraint_functions::AbstractVector{<:ConstraintFunction},
    save_folder::AbstractString,
    save_dd::Bool,
    p::ProgressMeter.Progress)

    # parallel evaluation of a generation
    ProgressMeter.next!(p)
    tmp = Distributed.pmap(
        x -> _optimization_engine(ini, act, actor_or_workflow, x, objective_functions, constraint_functions, save_folder, p.counter, save_dd),
        [X[k, :] for k in 1:size(X)[1]]
    )
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
