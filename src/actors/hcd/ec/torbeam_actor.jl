import TORBEAM
using Plots

#= =========== =#
#  ActorTORBEAM  #
#= =========== =#
@actor_parameters_struct ActorTORBEAM{T} begin
    backend::Switch{Symbol} = Switch{Symbol}([:fortran, :julia], "-", "TORBEAM implementation: the Fortran library (needs TORBEAM_DIR) or the pure-Julia one"; default=:fortran)
end

mutable struct ActorTORBEAM{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorTORBEAM{P}}
    torbeam_params::TORBEAM.TorbeamParams
end

function ActorTORBEAM(dd::IMAS.DD, par::FUSEparameters__ActorTORBEAM; kw...)
    logging_actor_init(ActorTORBEAM)
    par = OverrideParameters(par; kw...)
    # TorbeamParams gained `backend` in TORBEAM 1.1; keep working with older releases
    if hasfield(TORBEAM.TorbeamParams, :backend)
        torbeam_params = TORBEAM.TorbeamParams(; backend=par.backend)
    elseif par.backend == :fortran
        torbeam_params = TORBEAM.TorbeamParams()
    else
        error("ActorTORBEAM backend=:$(par.backend) needs TORBEAM.jl >= 1.1 (installed: $(pkgversion(TORBEAM)))")
    end
    return ActorTORBEAM(dd, par, torbeam_params)
end

"""
    ActorTORBEAM(dd::IMAS.DD, act::ParametersAllActors; kw...)

Performs electron cyclotron (EC) heating and current drive calculations using the 
TORBEAM ray-tracing code. TORBEAM provides detailed 3D beam propagation modeling
for EC waves, including absorption and current drive efficiency calculations.

The actor interfaces with the external TORBEAM code to perform sophisticated 
beam physics calculations that account for relativistic effects, mode conversion,
and realistic wave-plasma interactions.

!!! note

    With `backend=:fortran` (default) this requires the external TORBEAM library
    (`TORBEAM_DIR`); `backend=:julia` runs TORBEAM.jl's pure-Julia implementation.
    Reads data from `dd.ec_launchers` and equilibrium data, stores results in
    appropriate IMAS data structures.
"""
function ActorTORBEAM(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorTORBEAM(dd, act.ActorTORBEAM; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorTORBEAM)
    dd = actor.dd
    TORBEAM.run_torbeam(dd, actor.torbeam_params)
    return actor
end
