#= ============= =#
#  ActorSawteethSource  #
#= ============= =#
@actor_parameters_struct ActorSawteethSource{T} begin
    flat_factor::Entry{Float64} = Entry{Float64}("-", "Degree of flattening"; default=0.5, check=x -> @assert 0.0 <= x <= 1.0 "must be: 0.0 <= flat_factor <= 1.0")
    period::Entry{Float64} = Entry{Float64}("s", "Time after last sawteeth event to keep full flattening"; default=Inf, check=x -> @assert x >= 0.0 "period must be >= 0.0")
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

    Optional period keeps the flattening strength fully active for a
    specified time after the last sawteeth event and then ramps down linearly
    over an equal time interval.
"""
function ActorSawteethSource(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSawteethSource(dd, act.ActorSawteethSource; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _last_sawteeth_event(dd::IMAS.dd)
    diag = dd.sawteeth.diagnostics
    if ismissing(dd.sawteeth, :time) || ismissing(diag, :rho_tor_norm_inversion) || isempty(diag.rho_tor_norm_inversion)
        return nothing, nothing
    end
    @assert length(diag.rho_tor_norm_inversion) == length(dd.sawteeth.time) "sawteeth diagnostics length mismatch"

    i_last = IMAS.nearest_causal_time(dd.sawteeth.time, dd.global_time; bounds_error=false).index
    t_last = dd.sawteeth.time[i_last]
    rho_last = diag.rho_tor_norm_inversion[i_last]

    return t_last, rho_last
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
    if isinf(par.period)
        IMAS.sawteeth_source!(dd; par.flat_factor)
        return actor
    end

    t_last, rho_last = _last_sawteeth_event(dd)

    if isnothing(t_last)
        @info "ActorSawteethSource: `period` is finite but no previous sawteeth events found; no sawteeth source applied."
        IMAS.sawteeth_source!(dd, 0.0; par.flat_factor)
        return actor
    end

    Δt_since = max(0.0, dd.global_time - t_last)
    if Δt_since <= par.period
        activity = 1.0
    elseif Δt_since < 2.0 * par.period
        activity = 1.0 - (Δt_since - par.period) / par.period
    else
        activity = 0.0
    end

    IMAS.sawteeth_source!(dd, rho_last; flat_factor=par.flat_factor * activity)

    return actor
end
