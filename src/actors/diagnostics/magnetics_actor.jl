#= ============== =#
#  ActorMagnetics  #
#= ============== =#

@actor_parameters_struct ActorMagnetics{T} begin
    #== actor parameters ==#
    #== display and debugging parameters ==#
end

mutable struct ActorMagnetics{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorMagnetics{P}}
end

"""
    ActorMagnetics(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates magnetic field diagnostics from equilibrium data.

This actor computes the magnetic field measurements that would be observed by
magnetic diagnostics (flux loops, magnetic probes, etc.) based on the plasma
equilibrium. It populates the `dd.magnetics` IDS with synthetic diagnostic signals.

The actor calls `IMAS.magnetics!` to:
- Calculate magnetic flux measurements at diagnostic locations
- Compute magnetic field components for probes
- Generate time-dependent synthetic signals consistent with the equilibrium evolution

Key outputs:
- Magnetic flux loop signals in `dd.magnetics.flux_loop`
- Magnetic probe measurements in `dd.magnetics.bpol_probe`
- Time-dependent diagnostic data for validation and analysis

!!! note

    Requires `dd.equilibrium` to be populated with equilibrium data.
    Results are stored in `dd.magnetics`.
"""
function ActorMagnetics(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorMagnetics(dd, act.ActorMagnetics; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorMagnetics(dd::IMAS.dd, par::FUSEparameters__ActorMagnetics; kw...)
    logging_actor_init(ActorMagnetics)
    par = OverrideParameters(par; kw...)
    return ActorMagnetics(dd, par)
end

function _step(actor::ActorMagnetics)
    dd = actor.dd
    par = actor.par

    # Calculate magnetics diagnostics from equilibrium
    IMAS.magnetics!(dd.magnetics, dd.equilibrium)

    return actor
end
