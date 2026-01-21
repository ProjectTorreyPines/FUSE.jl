#= ============= =#
#  ActorSawteethSource  #
#= ============= =#
@actor_parameters_struct ActorSawteethSource{T} begin
    flat_factor::Entry{Float64} = Entry{Float64}("-", "Degree of flattening"; default=0.5, check=x -> @assert 0.0 <= x <= 1.0 "must be: 0.0 <= flat_factor <= 1.0")
    period::Entry{Float64} = Entry{Float64}("s", "Time after last sawteeth event to keep full flattening"; default=Inf, check=x -> @assert x >= 0.0 "period must be >= 0.0")
    qmin_desired::Entry{Float64} = Entry{Float64}("-", "Fallback qmin_desired when sawteeth diagnostics are missing"; default=1.0, check=x -> @assert x >= 0.0 "qmin_desired >= 0.0")
    ignore_before_time::Entry{Union{Nothing,Float64}} = Entry{Union{Nothing,Float64}}("s", "Ignore sawteeth diagnostics at times <= this value (nothing = use all history)"; default=nothing)
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

This actor modifies core sources using `IMAS.sawteeth_source!` to account for the redistribution
of heat and particles due to sawtooth crashes within the q=1 region. The sawtooth model
flattens source profiles from the core to the q=1 surface, simulating the mixing that occurs
during sawtooth crashes in tokamak plasmas.

# Parameters
- `flat_factor`: Degree of source flattening (0.0 = no flattening, 1.0 = full flattening)
- `period`: Controls time-dependent flattening behavior:
  - `Inf` (default): Always applies flattening based on current q-profile
  - Finite value: Uses sawteeth events recorded in `dd.sawteeth.diagnostics.rho_tor_norm_inversion`
    (typically written by ActorQED). Flattening is full strength for `period` seconds after
    the last event, then ramps linearly to zero over the next `period` seconds.
- `qmin_desired`: Fallback threshold used when no non-zero sawteeth events are found
- `ignore_before_time`: Ignore sawteeth diagnostics at times <= this value. If `nothing` (default),
    considers entire history. Only non-zero `rho_tor_norm_inversion` values are considered as sawteeth events.
"""
function ActorSawteethSource(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSawteethSource(dd, act.ActorSawteethSource; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _last_sawteeth_event(dd::IMAS.dd, ignore_before_time::Union{Nothing,Float64})
    diag = dd.sawteeth.diagnostics
    if ismissing(dd.sawteeth, :time) || ismissing(diag, :rho_tor_norm_inversion) || isempty(diag.rho_tor_norm_inversion)
        return nothing, nothing
    end
    @assert length(diag.rho_tor_norm_inversion) == length(dd.sawteeth.time) "sawteeth diagnostics length mismatch"
    info = IMAS.nearest_causal_time(dd.sawteeth.time, dd.global_time; bounds_error=false)
    if info.out_of_bounds
        return nothing, nothing
    end
    i_current = info.index

    # Find last non-zero inversion radius, filtered by time if ignore_before_time is set
    i_last = findlast(eachindex(diag.rho_tor_norm_inversion)) do i
        i <= i_current &&
        (ignore_before_time === nothing || dd.sawteeth.time[i] > ignore_before_time) &&
        diag.rho_tor_norm_inversion[i] > 0.0
    end

    if i_last === nothing
        return nothing, nothing
    end

    t_last = dd.sawteeth.time[i_last]
    rho_last = diag.rho_tor_norm_inversion[i_last]

    return t_last, rho_last
end

"""
    _step(actor::ActorSawteethSource)

Applies sawtooth source redistribution.

When `period = Inf`: Finds the last non-zero `rho_tor_norm_inversion` in the diagnostics
(optionally filtered by `ignore_before_time`) and applies flattening at that radius.
If no non-zero events are found, falls back to `qmin_desired` threshold.

When `period` is finite: Retrieves the most recent non-zero sawteeth event from
`dd.sawteeth` and applies time-decaying flattening. The flattening strength is:
- Full (`flat_factor`) for time ≤ `period` after the event
- Linearly ramping to zero for time between `period` and `2*period`
- Zero for time > `2*period` after the event
"""
function _step(actor::ActorSawteethSource)
    dd = actor.dd
    par = actor.par
    if isinf(par.period)
        _, rho_last = _last_sawteeth_event(dd, par.ignore_before_time)
        if rho_last === nothing
            # No non-zero sawteeth events found — use fallback qmin_desired
            IMAS.sawteeth_source!(dd; qmin_desired=par.qmin_desired, flat_factor=par.flat_factor)
        else
            # Use the last non-zero inversion radius
            IMAS.sawteeth_source!(dd, rho_last; par.flat_factor)
        end
        return actor
    end

    t_last, rho_last = _last_sawteeth_event(dd, par.ignore_before_time)

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
