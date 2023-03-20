#= == =#
#  IC  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorICsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Width of the deposition profile"; default=0.1)
    rho_0::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Radial location of the deposition profile"; default=0.0)
    current_efficiency::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "A/W", "Current drive efficiency"; default=0.125)
end

mutable struct ActorICsimple <: HCDAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorICsimple

    function ActorICsimple(dd::IMAS.dd, par::FUSEparameters__ActorICsimple; kw...)
        logging_actor_init(ActorICsimple)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorICsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the ion-cyclotron electron/ion energy deposition and current drive as a gaussian.

!!! note
    Reads data in `dd.ic_antennas` and stores data in `dd.core_sources`
"""
function ActorICsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorICsimple(kw...)
    actor = ActorICsimple(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorICsimple)
    dd = actor.dd
    par = actor.par

    n_antennas = length(dd.ic_antennas.antenna)
    _, width, rho_0, current_efficiency = same_length_vectors(1:n_antennas, par.width, par.rho_0, par.current_efficiency)

    for (idx, ica) in enumerate(dd.ic_antennas.antenna)
        eqt = dd.equilibrium.time_slice[]
        cp1d = dd.core_profiles.profiles_1d[]
        cs = dd.core_sources

        power_launched = @ddtime(ica.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
        area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

        # for FPP cases 80% to ions is reasonable (especially using minority heating)
        ion_electron_fraction_cp = fill(0.8, length(rho_cp))

        ne_vol = integrate(volume_cp, cp1d.electrons.density) / volume_cp[end]
        j_parallel = current_efficiency[idx] / eqt.boundary.geometric_axis.r / (ne_vol / 1e19) * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        source_index = IMAS.name_2_index(cs.source)[:ic]
        source = resize!(cs.source, "identifier.index" => source_index; allow_multiple_matches=true)
        gaussian_source(
            source,
            ica.name,
            source_index,
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