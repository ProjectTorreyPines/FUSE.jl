#= == =#
#  EC  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorECsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Width of the deposition profile"; default=0.1)
    rho_0::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Radial location of the deposition profile"; default=0.0)
    current_efficiency::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "A/W", "Current drive efficiency"; default=0.2)
end

mutable struct ActorECsimple <: HCDAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorECsimple
    width::AbstractVector{<:Real}
    rho_0::AbstractVector{<:Real}
    current_efficiency::AbstractVector{<:Real}
end

"""
    ActorECsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the EC electron energy deposition and current drive as a gaussian.

!!! note 
    Stores data in `dd.ec_launchers, dd.core_sources`
"""
function ActorECsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorECsimple(kw...)
    actor = ActorECsimple(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorECsimple(dd::IMAS.dd, par::FUSEparameters__ActorECsimple; kw...)
    logging_actor_init(ActorECsimple)
    par = par(kw...)
    n_launchers = length(dd.ec_launchers.beam)
    _, width, rho_0, current_efficiency = same_length_vectors(1:n_launchers, par.width, par.rho_0, par.current_efficiency)
    return ActorECsimple(dd, par, width, rho_0, current_efficiency)
end

function _step(actor::ActorECsimple)
    for (idx, ecl) in enumerate(actor.dd.ec_launchers.beam)
        eqt = actor.dd.equilibrium.time_slice[]
        cp1d = actor.dd.core_profiles.profiles_1d[]
        cs = actor.dd.core_sources

        power_launched = @ddtime(ecl.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
        area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

        ion_electron_fraction = 0.0

        ne_vol = integrate(volume_cp, cp1d.electrons.density) / volume_cp[end]
        j_parallel = actor.current_efficiency[idx] / eqt.boundary.geometric_axis.r / (ne_vol / 1e19) * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        source_index = IMAS.name_2_index(cs.source)[:ec]
        source = resize!(cs.source, "identifier.index" => source_index)
        gaussian_source(
            source,
            ecl.name,
            source_index,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction,
            actor.rho_0[idx],
            actor.width[idx],
            1.0;
            j_parallel=j_parallel
        )
    end
    return actor
end