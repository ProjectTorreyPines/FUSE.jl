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
    width::AbstractVector{<:Real}
    rho_0::AbstractVector{<:Real}
    current_efficiency::AbstractVector{<:Real}
end

"""
    ActorICsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the ion-cyclotron electron/ion energy deposition and current drive as a gaussian.

!!! note 
    Stores data in `dd.ic_antennas, dd.core_sources`
"""
function ActorICsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorICsimple(kw...)
    actor = ActorICsimple(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorICsimple(dd::IMAS.dd, par::FUSEparameters__ActorICsimple; kw...)
    logging_actor_init(ActorICsimple)
    par = par(kw...)
    n_antennas = length(dd.ic_antennas.antenna)
    _, width, rho_0, current_efficiency = same_length_vectors(1:n_antennas, par.width, par.rho_0, par.current_efficiency)
    return ActorICsimple(dd, par, width, rho_0, current_efficiency)
end

function _step(actor::ActorICsimple)
    for (idx, ica) in enumerate(actor.dd.ic_antennas.antenna)
        eqt = actor.dd.equilibrium.time_slice[]
        cp1d = actor.dd.core_profiles.profiles_1d[]
        cs = actor.dd.core_sources

        power_launched = @ddtime(ica.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
        area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

        ion_electron_fraction = 0.8 # for FPP cases 80% to ions is reasonable (especially using minority heating)

        ne_vol = integrate(volume_cp, cp1d.electrons.density) / volume_cp[end]
        j_parallel = actor.current_efficiency[idx] / eqt.boundary.geometric_axis.r / (ne_vol / 1e19) * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        isource = resize!(cs.source, "identifier.name" => ica.name)
        gaussian_source(
            isource,
            ica.name,
            5,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction,
            actor.rho_0[idx],
            actor.width[idx],
            1;
            j_parallel=j_parallel
        )
    end
    return actor
end