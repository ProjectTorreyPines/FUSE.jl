#= ======== =#
#  ActorHCD  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorHCD{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    ec_model::Switch{Symbol} = Switch{Symbol}([:ECsimple, :TORBEAM, :replay, :none], "-", "EC source actor to run"; default=:ECsimple)
    ic_model::Switch{Symbol} = Switch{Symbol}([:ICsimple, :replay, :none], "-", "IC source actor to run"; default=:ICsimple)
    lh_model::Switch{Symbol} = Switch{Symbol}([:LHsimple, :replay, :none], "-", "LH source actor to run"; default=:LHsimple)
    nb_model::Switch{Symbol} = Switch{Symbol}([:NBsimple, :RABBIT, :replay, :none], "-", "NB source actor to run"; default=:NBsimple)
    pellet_model::Switch{Symbol} = Switch{Symbol}([:Pelletsimple, :replay, :none], "-", "Pellet source actor to run"; default=:Pelletsimple)
    neutral_model::Switch{Symbol} = Switch{Symbol}([:neucg,:replay, :none], "-", "Pellet source actor to run"; default=:neucg)
end

mutable struct ActorHCD{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorHCD{P}
    act::ParametersAllActors{P}
    ec_actor::Union{ActorSimpleEC{D,P},ActorTORBEAM{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
    ic_actor::Union{ActorSimpleIC{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
    lh_actor::Union{ActorSimpleLH{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
    nb_actor::Union{ActorSimpleNB{D,P},ActorRABBIT{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
    #nb_actor::Union{ActorSimpleNB{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}

    pellet_actor::Union{ActorSimplePL{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
    neutral_actor::Union{ActorNeutralFueling{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
end

"""
    ActorHCD(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run HCD actors
"""
function ActorHCD(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorHCD(dd, act.ActorHCD, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorHCD(dd::IMAS.dd, par::FUSEparameters__ActorHCD, act::ParametersAllActors; kw...)
    logging_actor_init(ActorHCD)
    par = par(kw...)

    noop = ActorNoOperation(dd, act.ActorNoOperation)
    actor = ActorHCD(dd, par, act, noop, noop, noop, noop, noop, noop)

    @assert length(dd.pulse_schedule.ec.beam) == length(dd.ec_launchers.beam) "length(dd.pulse_schedule.ec.beam)=$(length(dd.pulse_schedule.ec.beam)) VS length(dd.ec_launchers.beam)=$(length(dd.ec_launchers.beam))"
    # fill missing EC launcher hardware details
    if par.ec_model in (:ECsimple, :TORBEAM)
        eqt = dd.equilibrium.time_slice[]
        for (k, ecb) in enumerate(dd.ec_launchers.beam)
            setup(ecb, eqt, dd.wall, act.ActorSimpleEC.actuator[k])
        end
    end
    if par.ec_model == :ECsimple
        actor.ec_actor = ActorSimpleEC(dd, act.ActorSimpleEC)
    elseif par.ec_model == :TORBEAM
        actor.ec_actor = ActorTORBEAM(dd, act.ActorTORBEAM)
    elseif par.ec_model == :replay
        actor.ec_actor = ActorReplay(dd, act.ActorReplay, actor.ec_actor)
    end

    @assert length(dd.pulse_schedule.ic.antenna) == length(dd.ic_antennas.antenna) "length(dd.pulse_schedule.ic.antenna)=$(length(dd.pulse_schedule.ic.antenna)) VS length(dd.ic_antennas.antenna)=$(length(dd.ic_antennas.antenna))"
    if par.ic_model == :ICsimple
        actor.ic_actor = ActorSimpleIC(dd, act.ActorSimpleIC)
    elseif par.ic_model == :replay
        actor.ic_actor = ActorReplay(dd, act.ActorReplay, actor.ic_actor)
    end

    @assert length(dd.pulse_schedule.lh.antenna) == length(dd.lh_antennas.antenna) "length(dd.pulse_schedule.lh.antenna)=$(length(dd.pulse_schedule.lh.antenna)) VS length(dd.lh_antennas.antenna)=$(length(dd.lh_antennas.antenna))"
    if par.lh_model == :LHsimple
        actor.lh_actor = ActorSimpleLH(dd, act.ActorSimpleLH)
    elseif par.lh_model == :replay
        actor.lh_actor = ActorReplay(dd, act.ActorReplay, actor.lh_actor)
    end

    @assert length(dd.pulse_schedule.nbi.unit) == length(dd.nbi.unit) "length(dd.pulse_schedule.nbi.unit)=$(length(dd.pulse_schedule.nbi.unit)) VS length(dd.nbi.unit)=$(length(dd.nbi.unit))"
    if par.nb_model == :NBsimple
        actor.nb_actor = ActorSimpleNB(dd, act.ActorSimpleNB)
    elseif par.nb_model == :RABBIT
        actor.nb_actor = ActorRABBIT(dd, act.ActorRABBIT)
    elseif par.nb_model == :replay
        actor.nb_actor = ActorReplay(dd, act.ActorReplay, actor.nb_actor)
    end

    @assert length(dd.pulse_schedule.pellet.launcher) == length(dd.pellets.launcher) "length(dd.pulse_schedule.pellet.launcher)=$(length(dd.pulse_schedule.pellet.launcher)) VS length(dd.pellets.launcher)=$(length(dd.pellets.launcher))"
    if par.pellet_model == :Pelletsimple
        actor.pellet_actor = ActorSimplePL(dd, act.ActorSimplePL)
    elseif par.pellet_model == :replay
        actor.pellet_actor = ActorReplay(dd, act.ActorReplay, actor.pellet_actor)
    end
    if par.neutral_model == :neucg
        actor.neutral_actor = ActorNeutralFueling(dd, act.ActorNeutralFueling)
    elseif par.neutral_model == :replay
        actor.neutral_actor = ActorReplay(dd, act.ActorReplay, actor.neutral_actor)
    end
    return actor

end

"""
    _step(actor::ActorHCD)

Runs through the selected HCD actor's step
"""
function _step(actor::ActorHCD)
    dd = actor.dd

    # Call IMAS.sources!(dd) since the most would expect sources to be consistent when coming out of this actor
    IMAS.sources!(dd)

    step(actor.ec_actor)
    step(actor.ic_actor)
    step(actor.lh_actor)
    step(actor.nb_actor)
    step(actor.pellet_actor)
    step(actor.neutral_actor)

    return actor
end

"""
    _finalize(actor::ActorHCD)

Finalizes the selected CHD actor's finalize
"""
function _finalize(actor::ActorHCD)

    finalize(actor.ec_actor)
    finalize(actor.ic_actor)
    finalize(actor.lh_actor)
    finalize(actor.nb_actor)
    finalize(actor.pellet_actor)
    finalize(actor.neutral_actor)
    return actor
end

function setup(ecb::IMAS.ec_launchers__beam, eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall, par::_FUSEparameters__ActorSimpleECactuator)
    # Estimate operating frequency and mode
    if ismissing(ecb.frequency, :data)
        resonance = IMAS.ech_resonance(eqt)
        ecb.frequency.time = [-Inf]
        ecb.frequency.data = [resonance.frequency]
        ecb.mode = resonance.mode == "X" ? -1 : 1
    end
    # Pick a reasonable launch location
    if ismissing(ecb.launching_position, :r) || ismissing(ecb.launching_position, :z)
        fw = IMAS.first_wall(wall)
        index = argmax(fw.r .+ fw.z)
        @ddtime(ecb.launching_position.r = fw.r[index])
        @ddtime(ecb.launching_position.z = fw.z[index])
    end
    if ismissing(ecb.launching_position, :phi)
        @ddtime(ecb.launching_position.phi = 0.0)
    end
    # beam properties
    if ismissing(ecb.phase, :angle) || ismissing(ecb.phase, :curvature)
        @ddtime(ecb.phase.angle = 0.0)
        @ddtime(ecb.phase.curvature = [0.0, 0.0])
    end
    if ismissing(ecb.spot, :angle) || ismissing(ecb.spot, :size)
        @ddtime(ecb.spot.angle = 0.0)
        @ddtime(ecb.spot.size = [0.0172, 0.0172])
    end
    # aiming based on rho0
    if (ismissing(ecb, :steering_angle_tor) || ismissing(ecb, :steering_angle_pol)) && !ismissing(par, :rho_0)
        launch_r = @ddtime(ecb.launching_position.r)
        launch_z = @ddtime(ecb.launching_position.z)
        resonance_layer = IMAS.ech_resonance_layer(eqt, IMAS.frequency(ecb))
        _, _, RHO_interpolant = IMAS.ρ_interpolant(eqt)
        rho_resonance_layer = RHO_interpolant.(resonance_layer.r, resonance_layer.z)
        index = resonance_layer.z .> eqt.global_quantities.magnetic_axis.z
        sub_index = argmin(abs.(rho_resonance_layer[index] .- par.rho_0))
        @ddtime(ecb.steering_angle_tor = 0.0)
        @ddtime(ecb.steering_angle_pol = atan(resonance_layer.z[index][sub_index] - launch_z, resonance_layer.r[index][sub_index] - launch_r) + pi)
    end
end

function _step(replay_actor::ActorReplay, actor::ActorSimpleEC, replay_dd::IMAS.dd)
    IMAS.copy_timeslice!(actor.dd.ec_launcher, replay_dd.ec_launcher, actor.dd.global_time)
    copy_source_timeslice!(actor.dd, replay_dd, :ec)
    return replay_actor
end

function _step(replay_actor::ActorReplay, actor::ActorSimpleIC, replay_dd::IMAS.dd)
    IMAS.copy_timeslice!(actor.dd.ic_antenna, replay_dd.ic_antenna, actor.dd.global_time)
    copy_source_timeslice!(actor.dd, replay_dd, :ic)
    return replay_actor
end

function _step(replay_actor::ActorReplay, actor::ActorSimpleLH, replay_dd::IMAS.dd)
    IMAS.copy_timeslice!(actor.dd.lh_antenna, replay_dd.lh_antenna, actor.dd.global_time)
    copy_source_timeslice!(actor.dd, replay_dd, :lh)
    return replay_actor
end

function _step(replay_actor::ActorReplay, actor::ActorSimpleNB, replay_dd::IMAS.dd)
    IMAS.copy_timeslice!(actor.dd.nbi, replay_dd.nbi, actor.dd.global_time)
    copy_source_timeslice!(actor.dd, replay_dd, :nbi)
    return replay_actor
end

function _step(replay_actor::ActorReplay, actor::ActorSimplePL, replay_dd::IMAS.dd)
    IMAS.copy_timeslice!(actor.dd.pellet, replay_dd.pellet, actor.dd.global_time)
    copy_source_timeslice!(actor.dd, replay_dd, :pellet)
    return replay_actor
end

function copy_source_timeslice!(dd::IMAS.dd, replay_dd::IMAS.dd, identifier::Symbol)
    for replay_source in findall(identifier, replay_dd.core_sources.source)
        source = resize!(dd.core_sources.source, :identifier, "identifier.name" => replay_source.name; wipe=false)
        IMAS.copy_timeslice!(source, replay_source, actor.dd.global_time)
    end
end
