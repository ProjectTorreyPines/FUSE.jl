#= === =#
#  NBI  #
#= === =#
Base.@kwdef mutable struct FUSEparameters__ActorNBIsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Width of the deposition profile"; default=0.3)
    rho_0::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Radial location of the deposition profile"; default=0.0)
    current_efficiency::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "A/W", "Current drive efficiency"; default=0.3)
end

mutable struct ActorNBIsimple <: HCDAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorNBIsimple
    width::AbstractVector{<:Real}
    rho_0::AbstractVector{<:Real}
    current_efficiency::AbstractVector{<:Real}
end

"""
    ActorNBIsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the NBI ion/electron energy deposition, particle source, rotation and current drive source with a super-gaussian.

!!! note 
    Stores data in `dd.nbi, dd.core_sources`
"""
function ActorNBIsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorNBIsimple(kw...)
    actor = ActorNBIsimple(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorNBIsimple(dd::IMAS.dd, par::FUSEparameters__ActorNBIsimple; kw...)
    logging_actor_init(ActorNBIsimple)
    par = par(kw...)
    n_beams = length(dd.nbi.unit)
    _, width, rho_0, current_efficiency = same_length_vectors(1:n_beams, par.width, par.rho_0, par.current_efficiency)
    return ActorNBIsimple(dd, par, width, rho_0, current_efficiency)
end

function _step(actor::ActorNBIsimple)
    for (idx, nbu) in enumerate(actor.dd.nbi.unit)
        eqt = actor.dd.equilibrium.time_slice[]
        cp1d = actor.dd.core_profiles.profiles_1d[]
        cs = actor.dd.core_sources

        beam_energy = @ddtime (nbu.energy.data)
        beam_mass = nbu.species.a
        power_launched = @ddtime(nbu.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
        area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

        ion_electron_fraction = IMAS.sivukhin_fraction(cp1d, beam_energy, beam_mass)

        beam_particles = power_launched / (beam_energy * constants.e)
        momentum_source =
            sin(nbu.beamlets_group[1].angle) * beam_particles * sqrt(2 * beam_energy * constants.e / beam_mass / constants.m_u) * beam_mass * constants.m_u

        ne_vol = integrate(volume_cp, cp1d.electrons.density) / volume_cp[end]
        j_parallel = actor.current_efficiency[idx] / eqt.boundary.geometric_axis.r / (ne_vol / 1e19) * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        isource = resize!(cs.source, "identifier.name" => nbu.name)
        gaussian_source(
            isource,
            nbu.name,
            2,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction,
            actor.rho_0[idx],
            actor.width[idx],
            2.0;
            electrons_particles=beam_particles,
            momentum_tor=momentum_source,
            j_parallel=j_parallel
        )
    end
    return actor
end