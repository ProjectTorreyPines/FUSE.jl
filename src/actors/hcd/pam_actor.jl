import PAM

#= ======== =#
#  ActorPAM  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorPAM{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    drift_model::Switch{Symbol} = Switch{Symbol}([:HPI2, :Parks, :none], "-", "drift model"; default=:Parks)
    time_from::Switch{Symbol} = Switch{Symbol}([:pulse_schedule, :pellet_time, :none], "-", "initialize time for the pellet calculations"; default=:none)
    time_step::Entry{Float64} = Entry{Float64}("-", "Time step [s]"; default=0.00001)
    time_end::Entry{Float64} = Entry{Float64}("-", "Time end [s]"; default=0.0008)
    Bt_dependance::Entry{Bool} = Entry{Bool}("-", "Enable Bt dependance"; default=true)
    density_update::Entry{Bool} = Entry{Bool}("-", "Update plasma density"; default=false)
end

mutable struct ActorPAM{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorPAM{P}}
    outputs::Union{Nothing,PAM.Pellet}
    function ActorPAM(dd::IMAS.dd{D}, par::FUSEparameters__ActorPAM{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorPAM)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par, nothing)
    end
end


"""
    ActorPAM(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the Pellet particle direction, ablation rate, density source deposition

!!! note

    Reads data in `pellet`, `pulse_schedule`, `equilibrium`, and stores data in `pellet`
"""
function ActorPAM(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPAM(dd, act.ActorPAM; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorPAM)
    dd = actor.dd
    par = actor.par

    if par.time_from == :none
        t0 = 0.0
    else
        error("ActorPAM: time_from=$(par.time_from) is not implemented yet (only :none is supported)")
    end

    # Ensure the simulated flight is long enough for the pellet to cross the plasma:
    # the launch-to-magnetic-axis distance sets a machine-independent lower bound on t_finish
    # (PAM returns early once the pellet is fully ablated). par.time_end acts as a floor.
    pellet = dd.pellets.time_slice[].pellet[1]
    fp = pellet.path_geometry.first_point
    sp = pellet.path_geometry.second_point
    flight_distance = sqrt((fp.r - sp.r)^2 + (fp.z - sp.z)^2)
    t_finish = max(par.time_end, 2.0 * flight_distance / pellet.velocity_initial)

    inputs = (
        t_start=t0,
        t_finish=t_finish,
        time_step=par.time_step,
        drift_model=par.drift_model,
        Bt_dependance=par.Bt_dependance,
        update_plasma=par.density_update
    )
    actor.outputs = PAM.run_PAM(dd; inputs...)

    return actor
end
"""
    _finalize(actor::ActorPAM)

Deposits the pellet particle source in `dd.core_sources`.

Following `ActorSimplePL`, the source is the time-averaged deposition set by the pellet
injection frequency: PAM resolves a single pellet flight, we ignore those timing details
and apply the resulting deposition profile at the launcher repetition rate.
"""
function _finalize(actor::ActorPAM)
    dd = actor.dd
    cs = dd.core_sources
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    output = actor.outputs

    @assert !isempty(dd.pulse_schedule.pellet.launcher) "ActorPAM requires dd.pulse_schedule.pellet.launcher to set the injection frequency"

    rho_cp = cp1d.grid.rho_tor_norm
    volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
    area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

    # PAM returns the electron density [particles/m^3] deposited by a single pellet, summed
    # over the pellet flight. trace_surfaces (used internally by PAM) produces one flux
    # surface per eqt.profiles_1d.psi, so the surface dimension is already on the
    # equilibrium rho_tor_norm grid; interpolate it onto the core_profiles grid.
    density_increment = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, vec(sum(output.density_source; dims=1))).(rho_cp)

    # PAM currently resolves a single pellet (run_PAM uses dd.pellets.time_slice[].pellet[1])
    ps = dd.pulse_schedule.pellet.launcher[1]
    name = isempty(dd.pellets.launcher) ? "pellet" : dd.pellets.launcher[1].name

    # time-averaged particle source rate [particles/(m^3 s)] = per-pellet deposition × injection frequency
    frequency = @ddtime(ps.frequency.reference)
    electrons_particles = density_increment .* frequency

    source = resize!(cs.source, :pellet, "identifier.name" => name; wipe=false)
    IMAS.new_source(
        source,
        source.identifier.index,
        name,
        rho_cp,
        volume_cp,
        area_cp;
        electrons_particles
    )

    return actor
end


