#= ================ =#
#  ActorSOL  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorSOL{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:box, :replay, :none], "-", "SOL actor to run"; default=:box)
    do_debug::Entry{Bool} = Entry{Bool}("-","Flag for debugging"; default = false)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot = false)
end

mutable struct ActorSOL{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSOL{P}}
    act::ParametersAllActors{P}
    SOL_actor::Union{Nothing,ActorSOLBox{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
end

"""
    ActorSOL(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run different SOL actors
"""
function ActorSOL(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSOL(dd, act.ActorSOL, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorSOL(dd::IMAS.dd, par::FUSEparameters__ActorSOL, act::ParametersAllActors; kw...)
    logging_actor_init(ActorSOL)
    par = OverrideParameters(par; kw...)

    noop = ActorNoOperation(dd, act.ActorNoOperation)
    actor = ActorSOL(dd, par, act, noop)

    if par.model == :box
        actor.SOL_actor = ActorSOLBox(dd, act.ActorSOLBox)
    elseif par.model == :replay
        actor.SOL_actor = ActorReplay(dd, act.ActorReplay, actor)
    end

    return actor
end

"""
    _step(actor::ActorSOL)

Clears and initializes data in ? for SOL actors to run properly, then calls the `step()` function of the selected SOL actor
"""
function _step(actor::ActorSOL)
    dd = actor.dd
    par = actor.par

    if par.model !== :none
        # initialize dd.edge_profiles.profiles_1d for equilibrium actors
        prepare(actor)
    end

    # step selected SOL actor
    step(actor.SOL_actor)

    return actor
end

"""
    _finalize(actor::ActorSOL)

Calls the `finalize()` function of the selected SOL actor and populates edge profile information
"""
function _finalize(actor::ActorSOL)
    dd = actor.dd
    par = actor.par

    # finalize selected SOL actor
    finalize(actor.SOL_actor)

    return actor
end

"""
    prepare(actor::ActorSOL)

Prepare `dd.edge_profiles.profiles_1d` to run SOL actors

  - TBD
"""
function prepare(actor::ActorSOL)
    dd = actor.dd
    par = actor.par
    act = actor.act

    ps = dd.pulse_schedule
    pc = ps.position_control

    # TBD

    return dd
end