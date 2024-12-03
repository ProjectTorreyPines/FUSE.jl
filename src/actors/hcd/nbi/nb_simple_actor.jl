#= ========== =#
#  Simple NBI  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorSimpleNB{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    ηcd_scale::Entry{Union{T,Vector{T}}} = Entry{Union{T,Vector{T}}}("-", "Scaling factor for nominal current drive efficiency"; default=1.0)
end

mutable struct ActorSimpleNB{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorSimpleNB{P}
    function ActorSimpleNB(dd::IMAS.dd{D}, par::FUSEparameters__ActorSimpleNB{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSimpleNB)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSimpleNB(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the NBI ion/electron energy deposition, particle source, rotation and current drive source with a super-gaussian.

NOTE: Current drive efficiency from GASC, based on "G. Tonon 'Current Drive Efficiency Requirements for an Attractive Steady-State Reactor'"

!!! note

    Reads data in `dd.nbi`, `dd.pulse_schedule` and stores data in `dd.core_sources`
"""
function ActorSimpleNB(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSimpleNB(dd, act.ActorSimpleNB; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorSimpleNB)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    cs = dd.core_sources

    R0 = eqt.boundary.geometric_axis.r
    rho_cp = cp1d.grid.rho_tor_norm
    volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
    area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

    for (k, (ps, nbu)) in enumerate(zip(dd.pulse_schedule.nbi.unit, dd.nbi.unit))
        power_launched = @ddtime(ps.power.reference)
        rho_0 = @ddtime(ps.deposition_rho_tor_norm.reference)
        width = @ddtime(ps.deposition_rho_tor_norm_width.reference)
        beam_energy = @ddtime(nbu.energy.data)
        beam_mass = nbu.species.a

        @ddtime(nbu.power_launched.data = power_launched)

        ion_electron_fraction_cp = IMAS.sivukhin_fraction(cp1d, beam_energy, beam_mass)

        electrons_particles = power_launched / (beam_energy * IMAS.mks.e)
        momentum_tor =
            power_launched * sin(nbu.beamlets_group[1].angle) * electrons_particles * sqrt(2.0 * beam_energy * IMAS.mks.e / beam_mass / IMAS.mks.m_u) * beam_mass * IMAS.mks.m_u

        ne20 = IMAS.interp1d(rho_cp, cp1d.electrons.density).(rho_0) / 1E20
        TekeV = IMAS.interp1d(rho_cp, cp1d.electrons.temperature).(rho_0) / 1E3

        if typeof(par.ηcd_scale) <: Vector
            ηcd_scale = par.ηcd_scale[k]
        else
            ηcd_scale = par.ηcd_scale
        end

        eta = ηcd_scale * TekeV * 0.025
        j_parallel = eta / R0 / ne20 * power_launched
        j_parallel *= sign(eqt.global_quantities.ip) .* (1 .- ion_electron_fraction_cp)

        source = resize!(cs.source, :nbi, "identifier.name" => nbu.name; wipe=false)
        shaped_source(
            source,
            nbu.name,
            source.identifier.index,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction_cp,
            ρ -> gaus(ρ, rho_0, width, 2.0);
            electrons_particles,
            momentum_tor,
            j_parallel
        )

        # add nbi fast ion particles source
        ion = resize!(source.profiles_1d[].ion, 1)[1]
        IMAS.ion_element!(ion, 1, nbu.species.a; fast=true)
        ion.particles = source.profiles_1d[].electrons.particles
        ion.particles_inside = source.profiles_1d[].electrons.particles_inside
        ion.fast_particles_energy = beam_energy

    end
    return actor
end