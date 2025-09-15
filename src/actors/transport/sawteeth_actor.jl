#= ============= =#
#  ActorSawteeth  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorSawteeth{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
end

mutable struct ActorSawteeth{D,P} <: AbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSawteeth{P}}
    function ActorSawteeth(dd::IMAS.dd{D}, par::FUSEparameters__ActorSawteeth{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSawteeth)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSawteeth(dd::IMAS.dd, act::ParametersAllActors; kw...)

Applies sawteeth reconnection to core_sources and core_profiles when q < 1.0.

This actor finds the last radial position where q < 1.0 and applies
sawteeth source and profile modifications using IMAS.sawteeth_source!
and IMAS.sawteeth_profiles! functions.
"""
function ActorSawteeth(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSawteeth(dd, act.ActorSawteeth; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorSawteeth)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]

    # update core_sources and core_profiles with sawteeth

    i_qdes = findlast(cp1d.q .<  1.0)
    if i_qdes !== nothing
        rho_qdes = cp1d.grid.rho_tor_norm[i_qdes]
        IMAS.sawteeth_source!(dd, rho_qdes)
    end

    # # sawteeth_profiles! must have slightly higher q to avoid big jumps in bootstrap current
    # i_qdes = findlast(cp1d.q .<  1.1)
    # if i_qdes !== nothing
    #     rho_qdes = cp1d.grid.rho_tor_norm[i_qdes]
    #     IMAS.sawteeth_profiles!(cp1d, rho_qdes)
    # end

    return actor
end