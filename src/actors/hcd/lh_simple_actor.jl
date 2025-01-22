#= == =#
#  LH  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorSimpleLH{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    actuator::ParametersVector{_FUSEparameters__ActorSimple{T}} = ParametersVector{_FUSEparameters__ActorSimple{T}}()
end

function Base.resize!(par::FUSEparameters__ActorSimpleLH, n::Int)
    # default rho_0 and width
    resize!(par.actuator, n, 0.8, 0.05)
end

mutable struct ActorSimpleLH{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorSimpleLH{P}
    function ActorSimpleLH(dd::IMAS.dd{D}, par::FUSEparameters__ActorSimpleLH{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSimpleLH)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSimpleLH(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the Lower-hybrid electron energy deposition and current drive as a gaussian.

NOTE: Current drive efficiency from GASC, based on "G. Tonon 'Current Drive Efficiency Requirements for an Attractive Steady-State Reactor'"

!!! note

    Reads data in `dd.lh_antennas`, `dd.pulse_schedule` and stores data in `dd.core_sources`
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
        power_launched = @ddtime(ps.power.reference)
        rho_0 = par.actuator[k].rho_0
        width = par.actuator[k].width
        ηcd_scale = par.actuator[k].ηcd_scale

        @ddtime(lha.power_launched.data = power_launched)

        ion_electron_fraction_cp = zeros(length(rho_cp))

        ne20 = IMAS.interp1d(rho_cp, cp1d.electrons.density).(rho_0) / 1E20
        TekeV = IMAS.interp1d(rho_cp, cp1d.electrons.temperature).(rho_0) / 1E3
        zeff = IMAS.interp1d(rho_cp, cp1d.zeff).(rho_0)

        eta = ηcd_scale * TekeV * 0.037 * B0 / (5.0 + zeff) / ne20^0.33
        j_parallel = eta / R0 / ne20 * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        source = resize!(cs.source, :lh, "identifier.name" => lha.name; wipe=false)
        shaped_source(
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
    end
    return actor
end