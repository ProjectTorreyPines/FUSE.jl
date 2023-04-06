#= === =#
#  NBI  #
#= === =#
Base.@kwdef mutable struct FUSEparameters__ActorNBIsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Width of the deposition profile"; default=0.3)
    rho_0::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Radial location of the deposition profile"; default=0.0)
end

mutable struct ActorNBIsimple <: HCDAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorNBIsimple

    function ActorNBIsimple(dd::IMAS.dd, par::FUSEparameters__ActorNBIsimple; kw...)
        logging_actor_init(ActorNBIsimple)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorNBIsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the NBI ion/electron energy deposition, particle source, rotation and current drive source with a super-gaussian.

NOTE: Current drive efficiency from GASC, based on "G. Tonon 'Current Drive Efficiency Requirements for an Attractive Steady-State Reactor'"

!!! note
    Reads data in `dd.nbi` and stores data in `dd.core_sources`
"""
function ActorNBIsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorNBIsimple(kw...)
    actor = ActorNBIsimple(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorNBIsimple)
    dd = actor.dd
    par = actor.par

    n_beams = length(dd.nbi.unit)
    _, width, rho_0 = same_length_vectors(1:n_beams, par.width, par.rho_0)

    for (idx, nbu) in enumerate(dd.nbi.unit)
        eqt = dd.equilibrium.time_slice[]
        cp1d = dd.core_profiles.profiles_1d[]
        cs = dd.core_sources

        beam_energy = @ddtime (nbu.energy.data)
        beam_mass = nbu.species.a
        power_launched = @ddtime(nbu.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
        area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

        ion_electron_fraction_cp = IMAS.sivukhin_fraction(cp1d, beam_energy, beam_mass)

        beam_particles = power_launched / (beam_energy * constants.e)
        momentum_source =
            sin(nbu.beamlets_group[1].angle) * beam_particles * sqrt(2.0 * beam_energy * constants.e / beam_mass / constants.m_u) * beam_mass * constants.m_u

        R0 = eqt.boundary.geometric_axis.r
        ne20 = IMAS.interp1d(rho_cp, cp1d.electrons.density).(rho_0[idx]) / 1E20
        TekeV = IMAS.interp1d(rho_cp, cp1d.electrons.temperature).(rho_0[idx]) / 1E3

        eta = TekeV * 0.025
        j_parallel = eta / R0 / ne20 * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        source_index = IMAS.name_2_index(cs.source)[:nbi]
        source = resize!(cs.source, "identifier.index" => source_index; allow_multiple_matches=true)
        gaussian_source(
            source,
            nbu.name,
            source_index,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction_cp,
            rho_0[idx],
            width[idx],
            2.0;
            electrons_particles=beam_particles,
            momentum_tor=momentum_source,
            j_parallel=j_parallel
        )
    end
    return actor
end