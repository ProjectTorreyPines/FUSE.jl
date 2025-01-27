#= ========== =#
#  Simple NBI  #
#= ========== =#
Base.@kwdef mutable struct _FUSEparameters__ActorSimpleNBactuator{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    ηcd_scale::Entry{T} = Entry{T}("-", "Scaling factor for nominal current drive efficiency"; default=1.0)
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.0, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.3, check=x -> @assert x >= 0.0 "must be: width > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ActorSimpleNB{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    actuator::ParametersVector{_FUSEparameters__ActorSimpleNBactuator{T}} = ParametersVector{_FUSEparameters__ActorSimpleNBactuator{T}}()
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
        beam_energy = @ddtime(ps.energy.reference)
        rho_0 = par.actuator[k].rho_0
        width = par.actuator[k].width
        ηcd_scale = par.actuator[k].ηcd_scale

        @ddtime(nbu.power_launched.data = power_launched)
        @ddtime(nbu.energy.data = beam_energy)
        beam_mass = nbu.species.a
        ion_electron_fraction_cp = IMAS.sivukhin_fraction(cp1d, beam_energy, beam_mass)

        if beam_energy > 0.0
            particles_per_second = power_launched / (beam_energy * IMAS.mks.e) # [1/s]
        else
            particles_per_second = 0.0
        end
        velocity = sqrt(2.0 * beam_energy * IMAS.mks.e / (beam_mass * IMAS.mks.m_u)) # [m/s]
        momentum_tor = sin(nbu.beamlets_group[1].angle) * particles_per_second * velocity * beam_mass * IMAS.mks.m_u # [kg*m/s^2] = [N]

        ne20 = IMAS.interp1d(rho_cp, cp1d.electrons.density).(rho_0) / 1E20
        TekeV = IMAS.interp1d(rho_cp, cp1d.electrons.temperature).(rho_0) / 1E3

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
            ρ -> IMAS.gaus(ρ, rho_0, width, 2.0);
            electrons_particles=particles_per_second,
            momentum_tor,
            j_parallel
        )

        # add nbi fast ion particles source
        source1d = source.profiles_1d[]
        ion = resize!(source1d.ion, 1)[1]
        IMAS.ion_element!(ion, 1, nbu.species.a; fast=true)
        ion.particles = source1d.electrons.particles
        ion.particles_inside = source1d.electrons.particles_inside
        ion.fast_particles_energy = beam_energy

    end
    return actor
end