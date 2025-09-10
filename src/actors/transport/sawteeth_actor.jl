#= ============= =#
#  ActorSawteeth  #
#= ============= =#
@actor_parameters_struct ActorSawteeth{T} begin
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

Applies sawtooth reconnection effects to plasma sources when the safety factor q < 1.0.

This actor finds the outermost radial position where q < 1.0 (the q=1 surface) and applies
sawtooth source modifications using `IMAS.sawteeth_source!` to account for the redistribution
of heat and particles due to sawtooth crashes within the q=1 region.

The sawtooth model affects the core sources by redistributing energy and particles
from the core to the q=1 surface, simulating the flattening of profiles that occurs
during sawtooth crashes in tokamak plasmas.
"""
function ActorSawteeth(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSawteeth(dd, act.ActorSawteeth; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorSawteeth)

Identifies q=1 surface and applies sawtooth source redistribution.

Finds the outermost radial grid point where q < 1.0 and calls `IMAS.sawteeth_source!`
to redistribute heat and particle sources within the q=1 region, simulating the
profile flattening effects of sawtooth crashes.
"""
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