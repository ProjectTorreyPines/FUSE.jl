#= ========= =#
#  Simple LH  #
#= ========= =#
Base.@kwdef mutable struct _FUSEparameters__ActorSimpleLHactuator{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    ηcd_scale::Entry{T} = Entry{T}("-", "Scaling factor for nominal current drive efficiency"; default=1.0)
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.8, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.05, check=x -> @assert x >= 0.0 "must be: width > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ActorSimpleLH{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    actuator::ParametersVector{_FUSEparameters__ActorSimpleLHactuator{T}} = ParametersVector{_FUSEparameters__ActorSimpleLHactuator{T}}()
end

mutable struct ActorSimpleLH{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSimpleLH{P}}
    function ActorSimpleLH(dd::IMAS.dd{D}, par::FUSEparameters__ActorSimpleLH{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSimpleLH)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSimpleLH(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the Lower-hybrid electron energy deposition and current drive as a gaussian.

NOTE: Current drive efficiency from GASC, based on "G. Tonon 'Current Drive Efficiency Requirements for an Attractive Steady-State Reactor'"

!!! note

    Reads data in `dd.lh_antennas`, `dd.pulse_schedule` and stores data in `dd.waves` and `dd.core_sources`
"""
function ActorSimpleLH(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSimpleLH(dd, act.ActorSimpleLH; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorSimpleLH)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    cs = dd.core_sources

    R0 = eqt.boundary.geometric_axis.r
    B0 = abs(IMAS.B0_geo(eqt))
    rho_cp = cp1d.grid.rho_tor_norm
    volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
    area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

    for (k, (ps, lha)) in enumerate(zip(dd.pulse_schedule.lh.antenna, dd.lh_antennas.antenna))
        τ_th = 0.01 # what's a good averating time here?
        power_launched = max(0.0, IMAS.smooth_beam_power(dd.pulse_schedule.lh.time, ps.power.reference, dd.global_time, τ_th))
        rho_0 = par.actuator[k].rho_0
        width = par.actuator[k].width
        ηcd_scale = par.actuator[k].ηcd_scale

        coherent_wave = resize!(dd.waves.coherent_wave, "identifier.antenna_name" => lha.name; wipe=false)

        @ddtime(lha.power_launched.data = power_launched)

        ion_electron_fraction_cp = zeros(length(rho_cp))

        ne20 = IMAS.interp1d(rho_cp, cp1d.electrons.density).(rho_0) / 1E20
        TekeV = IMAS.interp1d(rho_cp, cp1d.electrons.temperature).(rho_0) / 1E3
        zeff = IMAS.interp1d(rho_cp, cp1d.zeff).(rho_0)

        eta = ηcd_scale * TekeV * 0.037 * B0 / (5.0 + zeff) / ne20^0.33
        j_parallel = eta / R0 / ne20 * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        source = resize!(cs.source, :lh, "identifier.name" => lha.name; wipe=false)
        shaped_source!(
            source,
            lha.name,
            source.identifier.index,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction_cp,
            ρ -> IMAS.gaus(ρ, rho_0, width, 1.0);
            j_parallel
        )

        # populate waves IDS
        resize!(coherent_wave.profiles_1d)
        populate_wave1d_from_source1d!(coherent_wave.profiles_1d[], source.profiles_1d[])
    end

    return actor
end