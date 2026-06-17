#= ================== =#
#  ActorInterferometer  #
#= ================== =#

@actor_parameters_struct ActorInterferometer{T} begin
    #== actor parameters ==#
    n_points::Entry{Int} = Entry{Int}("-", "Number of integration points along line of sight"; default=100)
    #== display and debugging parameters ==#
end

mutable struct ActorInterferometer{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorInterferometer{P}}
end

"""
    ActorInterferometer(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates synthetic interferometer measurements from equilibrium and core profiles.

This actor computes line-averaged electron density measurements that would be observed by
interferometer diagnostics based on the plasma equilibrium and core profiles. It populates
the `dd.interferometer` IDS with synthetic diagnostic signals.

The actor calls `IMAS.interferometer!` to:
- Calculate line-integrated electron density along each interferometer line of sight
- Compute line-averaged density for each channel
- Generate time-dependent synthetic signals consistent with the equilibrium and profile evolution

Key physics:
- Integration of electron density along interferometer chords
- Proper handling of line-of-sight geometry (including multi-segment paths)
- Time-dependent measurements across all equilibrium time slices

Key outputs:
- Line-averaged electron density in `dd.interferometer.channel[:].n_e_line_average`
- Time-dependent diagnostic data for validation and comparison with experimental data

Parameters:
- `n_points`: Number of integration points along each line of sight (default: 100)

!!! note

    Requires `dd.equilibrium` and `dd.core_profiles` to be populated.
    The `dd.interferometer` must have channels defined with `line_of_sight` geometry.
    Results are stored in `dd.interferometer.channel[:].n_e_line_average`.
"""
function ActorInterferometer(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorInterferometer(dd, act.ActorInterferometer; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorInterferometer(dd::IMAS.dd, par::FUSEparameters__ActorInterferometer; kw...)
    logging_actor_init(ActorInterferometer)
    par = OverrideParameters(par; kw...)
    return ActorInterferometer(dd, par)
end

function _step(actor::ActorInterferometer)
    dd = actor.dd
    par = actor.par

    # Calculate interferometer diagnostics from equilibrium and core profiles
    IMAS.interferometer!(dd.interferometer, dd.equilibrium, dd.core_profiles; n_points=par.n_points)

    return actor
end
