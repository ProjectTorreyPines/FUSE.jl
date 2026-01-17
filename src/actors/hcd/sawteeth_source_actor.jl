#= ============= =#
#  ActorSawteethSource  #
#= ============= =#
@actor_parameters_struct ActorSawteethSource{T} begin
    flat_factor::Entry{Float64} = Entry{Float64}("-", "Degree of flattening"; default=0.5, check=x -> @assert 0.0 <= x <= 1.0 "must be: 0.0 <= flat_factor <= 1.0")
end

mutable struct ActorSawteethSource{D,P} <: AbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSawteethSource{P}}
    function ActorSawteethSource(dd::IMAS.dd{D}, par::FUSEparameters__ActorSawteethSource{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSawteethSource)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSawteethSource(dd::IMAS.dd, act::ParametersAllActors; kw...)

Applies sawtooth reconnection effects to plasma sources when the safety factor q < 1.0.

This actor finds the outermost radial position where q < 1.0 (the q=1 surface) and applies
sawtooth source modifications using `IMAS.sawteeth_source!` to account for the redistribution
of heat and particles due to sawtooth crashes within the q=1 region.

The sawtooth model affects the core sources by redistributing energy and particles
from the core to the q=1 surface, simulating the flattening of profiles that occurs
during sawtooth crashes in tokamak plasmas.
"""
function ActorSawteethSource(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSawteethSource(dd, act.ActorSawteethSource; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorSawteethSource)

Identifies q=1 surface and applies sawtooth source redistribution.

Finds the outermost radial grid point where q < 1.0 and calls `IMAS.sawteeth_source!`
to redistribute heat and particle sources within the q=1 region, simulating the
profile flattening effects of sawtooth crashes.
"""
function _step(actor::ActorSawteethSource)
    dd = actor.dd
    par = actor.par
    IMAS.sawteeth_source!(dd; par.flat_factor)
    return actor
end