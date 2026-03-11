import TORBEAM
using Plots

#= =========== =#
#  ActorTORBEAM  #
#= =========== =#
@actor_parameters_struct ActorTORBEAM{T} begin
end

mutable struct ActorTORBEAM{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorTORBEAM{P}}
    torbeam_params::TORBEAM.TorbeamParams
end

function ActorTORBEAM(dd::IMAS.dd, par::FUSEparameters__ActorTORBEAM; kw...)
    logging_actor_init(ActorTORBEAM)
    par = OverrideParameters(par; kw...)
    return ActorTORBEAM(dd, par, TORBEAM.TorbeamParams())
end

"""
    ActorTORBEAM(dd::IMAS.dd, act::ParametersAllActors; kw...)

Performs electron cyclotron (EC) heating and current drive calculations using the 
TORBEAM ray-tracing code. TORBEAM provides detailed 3D beam propagation modeling
for EC waves, including absorption and current drive efficiency calculations.

The actor interfaces with the external TORBEAM code to perform sophisticated 
beam physics calculations that account for relativistic effects, mode conversion,
and realistic wave-plasma interactions.

!!! note

    Requires TORBEAM external code. Reads data from `dd.ec_launchers` and 
    equilibrium data, stores results in appropriate IMAS data structures.
"""
function ActorTORBEAM(dd::IMAS.dd, act::ParametersAllActors; kw...)
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
