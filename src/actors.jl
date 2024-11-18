import AbstractTrees
import OrderedCollections

abstract type AbstractActor{D,P} end

abstract type CompoundAbstractActor{D,P} <: AbstractActor{D,P} end

abstract type SingleAbstractActor{D,P} <: AbstractActor{D,P} end

function logging_actor_init(typeof_actor::Type{<:AbstractActor}, args...; kw...)
    return logging(Logging.Debug, :actors, "$(name(typeof_actor)) @ init")
end

#= =========== =#
#  name & group #
#= =========== =#
function name(actor::AbstractActor; remove_Actor::Bool=true)
    return name(typeof(actor); remove_Actor)
end

function name(typeof_actor::Type{<:AbstractActor}; remove_Actor::Bool=true)
    if remove_Actor
        return string(split(replace(string(typeof_actor), r"^FUSE\.Actor" => ""), "{")[1])
    else
        return string(split(replace(string(typeof_actor), r"^FUSE\." => ""), "{")[1])
    end
end

"""
    group_name(typeof_actor::Type{<:AbstractActor})

Returns the group to which the actor belongs to.

NOTE: This is based on the name of the folder in which the .jl file that defines the actor is contained.
"""
function group_name(typeof_actor::Type{<:AbstractActor})
    dd = IMAS.dd()
    act = ParametersActors()
    actor_name = name(typeof_actor; remove_Actor=false)
    if typeof_actor <: CompoundAbstractActor
        which_output = string(@which typeof_actor(dd, getproperty(act, Symbol(actor_name)), act))
    elseif typeof_actor <: SingleAbstractActor
        which_output = string(@which typeof_actor(dd, getproperty(act, Symbol(actor_name))))
    end
    folder_name = split(split(which_output, "@")[end], "/")[end-1]
    return folder_name
end

#= =============== =#
#  switch_get_from  #
#= =============== =#
"""
    switch_get_from(quantity::Symbol)

Switch to pick form which IDS `quantity` comes from
"""
function switch_get_from(quantity::Symbol)::Switch{Symbol}
    txt = "Take $quantity from this IDS"
    if quantity == :ip
        swch = Switch{Symbol}([:core_profiles, :equilibrium, :pulse_schedule], "-", txt)
    elseif quantity == :vacuum_r0_b0
        swch = Switch{Symbol}([:equilibrium, :pulse_schedule], "-", txt; default=:pulse_schedule)
    elseif quantity == :vloop
        swch = Switch{Symbol}([:core_profiles, :equilibrium, :pulse_schedule, :controllers__ip], "-", txt)
    elseif quantity == :βn
        swch = Switch{Symbol}([:core_profiles, :equilibrium], "-", txt)
    elseif quantity == :ne_ped
        swch = Switch{Symbol}([:core_profiles, :summary, :pulse_schedule], "-", txt)
    elseif quantity == :zeff_ped
        swch = Switch{Symbol}([:core_profiles, :summary, :pulse_schedule], "-", txt)
    else
        error("`$quantity` not supported in switch_get_from()")
    end
    return swch
end

#= ==================== =#
#  ParametersActor time  #
#= ==================== =#
function SimulationParameters.global_time(par::ParametersActor)
    return getfield(par, :_time)
end

function SimulationParameters.global_time(par::ParametersActor, time::Float64)
    return setfield!(par, :_time, time)
end

function SimulationParameters.time_range(par::ParametersActor)
    return missing
end

#= ==== =#
#  step  #
#= ==== =#
"""
    step(actor::T, args...; kw...) where {T<:AbstractActor}

Calls `_step(actor)`

This is where the main part of the actor calculation gets done
"""
function step(actor::T, args...; kw...) where {T<:AbstractActor}
    global_time(actor.par, global_time(actor.dd)) # syncrhonize global_time of `par` with the one of `dd`
    if !actor_logging(actor.dd)
        _step(actor, args...; kw...)::T
    else
        timer_name = name(actor)
        TimerOutputs.reset_timer!(timer_name)
        TimerOutputs.@timeit timer timer_name begin
            memory_time_tag(actor, "step IN")
            logging(Logging.Info, :actors, " "^workflow_depth(actor.dd) * "$(name(actor))")
            enter_workflow(actor)
            try
                s_actor = _step(actor, args...; kw...)
                @assert s_actor === actor "`$(T)._step(actor)` should return the same actor that is input to the function"
            finally
                exit_workflow(actor)
            end
            memory_time_tag(actor, "step OUT")
        end
    end
    return actor
end

#= ======== =#
#  finalize  #
#= ======== =#
function _finalize(actor::AbstractActor)
    return actor
end

"""
    finalize(actor::T) where {T<:AbstractActor}

Calls `_finalize(actor)`

This is typically used to update `dd` to whatever the actor has calculated at the `step` function
"""
function finalize(actor::T)::T where {T<:AbstractActor}
    if !actor_logging(actor.dd)
        _finalize_and_freeze_onetime_expressions(actor)::T
    else
        timer_name = name(actor)
        TimerOutputs.@timeit timer timer_name begin
            memory_time_tag(actor, "finalize IN")
            logging(Logging.Debug, :actors, " "^workflow_depth(actor.dd) * "$(name(actor)) @finalize")
            f_actor = _finalize_and_freeze_onetime_expressions(actor)
            @assert f_actor === actor "`$(T)._finalize(actor)` should return the same actor that is input to the function"
            memory_time_tag(actor, "finalize OUT")
        end
    end
    return actor
end

function _finalize_and_freeze_onetime_expressions(actor::T) where {T<:AbstractActor}
    f_actor = _finalize(actor)
    # freeze onetime expressions (ie. grids)
    while !isempty(IMAS.expression_onetime_weakref)
        idsw = pop!(IMAS.expression_onetime_weakref, first(keys(IMAS.expression_onetime_weakref)))
        if idsw.value !== nothing
            # println("Freeze $(typeof(actor)): $(IMAS.location(idsw.value))")
            IMAS.freeze!(idsw.value)
        end
    end
    return f_actor
end

@recipe function plot_actor(actor::AbstractActor, args...)
    return error("No plot recipe defined for ator $(typeof(actor))")
end

#= ============= =#
#  actor_logging  #
#= ============= =#
function actor_logging(dd::IMAS.dd)
    aux = getfield(dd, :_aux)
    if :fuse_actor_logging ∉ keys(aux)
        aux[:fuse_actor_logging] = true
    end
    return aux[:fuse_actor_logging]
end

function actor_logging(dd::IMAS.dd, value::Bool)
    aux = getfield(dd, :_aux)
    old_value = actor_logging(dd)
    aux[:fuse_actor_logging] = value
    return old_value
end

#= ==== =#
#  show  #
#= ==== =#
function Base.show(io::IO, ::MIME"text/plain", actor::AbstractActor)
    actorname = replace(string(typeof(actor)), "{Float64, Float64}" => "", r"^FUSE." => "")
    fields = join([fieldname for fieldname in fieldnames(typeof(actor))], ", ")
    return print(io, "$actorname($fields)")
end

#= ======== =#
#  workflow  #
#= ======== =#
mutable struct Workflow
    _name::String
    _flow::OrderedCollections.OrderedDict{Tuple,Workflow}
end
Workflow(name::String) = Workflow(name, OrderedCollections.OrderedDict{String,Workflow}())

function AbstractTrees.printnode(io::IO, workflow::Workflow)
    return printstyled(io, workflow._name)
end

function AbstractTrees.children(workflow::Workflow)
    return values(workflow._flow)
end

function Base.show(io::IO, ::MIME"text/plain", workflow::Workflow; maxdepth::Int=1000, kwargs...)
    return AbstractTrees.print_tree(io, workflow; maxdepth, kwargs...)
end

function enter_workflow(actor::AbstractActor)
    aux = _aux_workflow(actor.dd)
    h = goto_worflow_depth(aux[:fuse_workflow], aux[:fuse_workflow_depth])
    aux[:fuse_workflow_count] += 1
    aux[:fuse_workflow_depth] += 1
    return h[(name(actor), aux[:fuse_workflow_count])] = Workflow(name(actor))
end

function exit_workflow(actor::AbstractActor)
    return _aux_workflow(actor.dd)[:fuse_workflow_depth] -= 1
end

function goto_worflow_depth(workflow::Workflow, depth::Int)
    h = workflow._flow
    for k in 1:depth
        h = collect(values(h))[end]._flow
    end
    return h
end

function workflow_depth(dd::IMAS.dd)
    return _aux_workflow(dd)[:fuse_workflow_depth]
end

function workflow(dd::IMAS.dd)
    return _aux_workflow(dd)[:fuse_workflow]
end

function _aux_workflow(dd::IMAS.dd)
    aux = getfield(dd, :_aux)
    if :fuse_workflow ∉ keys(aux)
        aux[:fuse_workflow] = Workflow("Main")
        aux[:fuse_workflow_depth] = 0
        aux[:fuse_workflow_count] = 0
    end
    return aux
end
