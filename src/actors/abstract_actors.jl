abstract type AbstractActor end
abstract type FacilityAbstractActor <: AbstractActor end
abstract type ReactorAbstractActor <: AbstractActor end
abstract type HCDAbstractActor <: AbstractActor end
abstract type PlasmaAbstractActor <: AbstractActor end

function logging_actor_init(typeof_actor::DataType, args...; kw...)
    logging(Logging.Debug, :actors, "$typeof_actor @ init")
end

#= ==== =#
#  step  #
#= ==== =#
"""
    step(actor::T, args...; kw...) where {T<:AbstractActor}

Calls `_step(actor)`

This is where the calculation is done.
"""
function step(actor::T, args...; kw...) where {T<:AbstractActor}
    logging(Logging.Info, :actors, "$(typeof(actor)) @ step")
    timer_name = replace(string(typeof(actor).name.name),r"^Actor" => "")
    TimerOutputs.reset_timer!(timer_name)
    TimerOutputs.@timeit timer timer_name begin
        s = _step(actor, args...; kw...)
        @assert s === actor "_step should return the same actor (check if it is actor at all)"
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

This is typically used to update `dd` to whatever the actor has calculated at the `step` function.
"""
function finalize(actor::T) where {T<:AbstractActor}
    logging(Logging.Debug, :actors, "$(typeof(actor)) @finalize")
    s = _finalize(actor)
    @assert s === actor "_finalize should return the same actor (check if it is actor at all)"
    
    # freeze onetime expressions (ie. grids)
    while !isempty(IMASDD.expression_onetime_weakref)
        idsw = pop!(IMASDD.expression_onetime_weakref)
        if idsw.value !== nothing
            # println("Freeze $(typeof(actor)): $(IMAS.location(idsw.value))")
            IMAS.freeze(idsw.value)
        end
    end
    
    return actor
end
