#= ============ =#
#  ActorSolovev  #
#= ============ =#
include("solovev_actor.jl")

#= =========== =#
#  ActorCHEASE  #
#= =========== =#
include("chease_actor.jl")

#= ================ =#
#  ActorEquilibrium  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorEquilibrium{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model = Switch(Symbol, [:Solovev, :CHEASE], "-", "Equilibrium actor to run"; default=:Solovev)
end

mutable struct ActorEquilibrium <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorEquilibrium
    eq_actor::Union{ActorSolovev,ActorCHEASE}
end

"""
    ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple equilibrium actors
"""
function ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorEquilibrium(kw...)
    actor = ActorEquilibrium(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEquilibrium(dd::IMAS.dd, par::FUSEparameters__ActorEquilibrium, act::ParametersAllActors; kw...)
    logging_actor_init(ActorEquilibrium)
    par = par(kw...)
    if par.model == :Solovev
        eq_actor = ActorSolovev(dd, act.ActorSolovev)
    elseif par.model == :CHEASE
        eq_actor = ActorCHEASE(dd, act.ActorCHEASE)
    else
        error("ActorEquilibrium: model = $(par.model) is unknown")
    end
    return ActorEquilibrium(dd, par, eq_actor)
end

"""
    prepare(dd::IMAS.dd, :ActorEquilibrium, act::ParametersAllActors; kw...)

Prepare dd to run ActorEquilibrium
* call prapare function of the different equilibrium models
"""
function prepare(dd::IMAS.dd, ::Type{Val{:ActorEquilibrium}}, act::ParametersAllActors; kw...)
    par = act.ActorEquilibrium(kw...)
    if par.model == :Solovev
        return prepare(dd, :ActorSolovev, act; kw...)
    elseif par.model == :CHEASE
        return prepare(dd, :ActorCHEASE, act; kw...)
    else
        error("ActorEquilibrium: model = $(par.model) is unknown")
    end
end

"""
    step(actor::ActorEquilibrium)

Runs through the selected equilibrium actor's step
"""
function _step(actor::ActorEquilibrium)
    step(actor.eq_actor)
end

"""
    finalize(actor::ActorEquilibrium)

Finalizes the selected equilibrium actor
"""
function _finalize(actor::ActorEquilibrium)
    finalize(actor.eq_actor)
end