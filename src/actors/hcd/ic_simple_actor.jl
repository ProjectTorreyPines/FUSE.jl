#= == =#
#  IC  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorICsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width::Entry{Union{T,AbstractVector{T}}} = Entry{Union{T,AbstractVector{T}}}("-", "Width of the deposition profile"; default=0.1)
    rho_0::Entry{Union{T,AbstractVector{T}}} = Entry{Union{T,AbstractVector{T}}}("-", "Radial location of the deposition profile"; default=0.0)
end

mutable struct ActorICsimple{D,P} <: HCDAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorICsimple{P}
    function ActorICsimple(dd::IMAS.dd{D}, par::FUSEparameters__ActorICsimple{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorICsimple)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorICsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the ion-cyclotron electron/ion energy deposition and current drive as a gaussian.

NOTE: Current drive efficiency from GASC, based on "G. Tonon 'Current Drive Efficiency Requirements for an Attractive Steady-State Reactor'"

!!! note
    Reads data in `dd.ic_antennas` and stores data in `dd.core_sources`
"""
function ActorICsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorICsimple(dd, act.ActorICsimple; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorICsimple)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    cs = dd.core_sources

    R0 = eqt.boundary.geometric_axis.r
    beta_tor = eqt.global_quantities.beta_tor
    rho_cp = cp1d.grid.rho_tor_norm
    volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
    area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

    n_antennas = length(dd.ic_antennas.antenna)
    _, width, rho_0 = same_length_vectors(1:n_antennas, par.width, par.rho_0)

    for (idx, ica) in enumerate(dd.ic_antennas.antenna)
        power_launched = @ddtime(dd.pulse_schedule.ic.antenna[idx].power.reference.data)
        @ddtime(ica.power_launched.data = power_launched)

        # for FPP cases 80% to ions is reasonable (especially using minority heating)
        ion_electron_fraction_cp = fill(0.8, length(rho_cp))

        ne20 = IMAS.interp1d(rho_cp, cp1d.electrons.density).(rho_0[idx]) / 1E20
        TekeV = IMAS.interp1d(rho_cp, cp1d.electrons.temperature).(rho_0[idx]) / 1E3
        zeff = IMAS.interp1d(rho_cp, cp1d.zeff).(rho_0[idx])

        eta = TekeV * 0.063 / (2.0 + zeff) / (1.0 + 0.5 * beta_tor)
        j_parallel = eta / R0 / ne20 * power_launched
        j_parallel *= sign(eqt.global_quantities.ip) .* ion_electron_fraction_cp

        source = resize!(cs.source, :ic, "identifier.name" => ica.name; wipe=false)
        gaussian_source(
            source,
            ica.name,
            source.identifier.index,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction_cp,
            rho_0[idx],
            width[idx],
            1.0;
            j_parallel=j_parallel
        )
    end
    return actor
end