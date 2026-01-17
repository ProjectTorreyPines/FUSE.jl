#= ============= =#
#  ActorSawteeth  #
#= ============= =#
@actor_parameters_struct ActorSawteeth{T} begin
    flat_factor::Entry{Float64} = Entry{Float64}("-", "Degree of flattening"; default=0.5, check=x -> @assert 0.0 <= x <= 1.0 "must be: 0.0 <= flat_factor <= 1.0")
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
    par = actor.par
    IMAS.sawteeth_source!(dd; par.flat_factor)
    return actor
end