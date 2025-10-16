#= ================== =#
#  ActorCoreTransport  #
#= ================== =#
@actor_parameters_struct ActorCoreTransport{T} begin
    model::Switch{Symbol} = Switch{Symbol}([:FluxMatcher, :EPEDProfiles, :replay, :none], "-", "Transport actor to run"; default=:FluxMatcher)
end

mutable struct ActorCoreTransport{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorCoreTransport{P}}
    act::ParametersAllActors{P}
    tr_actor::Union{ActorFluxMatcher{D,P},ActorEPEDprofiles{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
end

"""
    ActorCoreTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a unified interface to run core transport evolution models.

This compound actor manages different approaches to core transport evolution:

Transport model options:
- `:FluxMatcher`: Self-consistent flux-matching transport evolution using turbulent and 
  neoclassical models to evolve temperature and density profiles
- `:EPEDProfiles`: Use EPED model predictions for pedestal and core profiles
- `:replay`: Replay profiles from experimental data or previous simulations
- `:none`: No core transport evolution (fixed profiles)

The selected model determines how the core plasma profiles (temperature, density, rotation)
evolve in response to heating, particle sources, and transport processes.
"""
function ActorCoreTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCoreTransport(dd, act.ActorCoreTransport, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCoreTransport(dd::IMAS.dd, par::FUSEparameters__ActorCoreTransport, act::ParametersAllActors; kw...)
    logging_actor_init(ActorCoreTransport)
    par = OverrideParameters(par; kw...)

    noop = ActorNoOperation(dd, act.ActorNoOperation)
    actor = ActorCoreTransport(dd, par, act, noop)

    if par.model == :FluxMatcher
        actor.tr_actor = ActorFluxMatcher(dd, act.ActorFluxMatcher, act)
    elseif par.model == :EPEDProfiles
        actor.tr_actor = ActorEPEDprofiles(dd, act.ActorEPEDprofiles, act)
    elseif par.model == :replay
        actor.tr_actor = ActorReplay(dd, act.ActorReplay, actor)
    end

    return actor
end

"""
    step(actor::ActorCoreTransport)

Runs through the selected core transport actor step
"""
function _step(actor::ActorCoreTransport)
    step(actor.tr_actor)
    return actor
end

"""
    _finalize(actor::ActorCoreTransport)

Finalizes the selected core transport actor finalize
"""
function _finalize(actor::ActorCoreTransport)
    finalize(actor.tr_actor)
    return actor
end

function _step(replay_actor::ActorReplay, actor::ActorCoreTransport, replay_dd::IMAS.dd)
    dd = actor.dd

    time0 = dd.global_time
    cp1d = dd.core_profiles.profiles_1d[time0]
    replay_cp1d = replay_dd.core_profiles.profiles_1d[time0]
    rho = cp1d.grid.rho_tor_norm

    # here we purposely set rho_nml == rho_ped
    # Here blend_core_edge() will connect the replay core to the edge,
    # but will allow for a discontinuity of the profiles.
    # It is the role of the ActorPedestal to do the actual blending between rho_nml and rho_ped
    rho_nml = actor.act.ActorFluxMatcher.rho_transport[end]
    rho_ped = actor.act.ActorFluxMatcher.rho_transport[end]

    # densities
    cp1d.electrons.density_thermal = IMAS.blend_core_edge(replay_cp1d.electrons.density_thermal, cp1d.electrons.density_thermal, rho, rho_nml, rho_ped; method=:scale)
    for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
        if !ismissing(ion, :density_thermal)
            ion.density_thermal = IMAS.blend_core_edge(replay_ion.density_thermal, ion.density_thermal, rho, rho_nml, rho_ped; method=:scale)
        end
    end
    IMAS.scale_ion_densities_to_target_zeff!(cp1d, replay_cp1d.zeff)

    # temperatures
    cp1d.electrons.temperature = IMAS.blend_core_edge(replay_cp1d.electrons.temperature, cp1d.electrons.temperature, rho, rho_nml, rho_ped)
    for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
        if !ismissing(ion, :temperature)
            ion.temperature = IMAS.blend_core_edge(replay_ion.temperature, ion.temperature, rho, rho_nml, rho_ped)
        end
    end

    # rotation
    cp1d.rotation_frequency_tor_sonic = replay_cp1d.rotation_frequency_tor_sonic

    return replay_actor
end