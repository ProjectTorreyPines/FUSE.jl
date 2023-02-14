import NumericalIntegration: integrate

#= === =#
#  NBI  #
#= === =#
Base.@kwdef mutable struct FUSEparameters__ActorNBIsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width = Entry(Union{Real,AbstractVector{<:Real}}, "-", "Width of the deposition profile"; default=0.3)
    rho_0 = Entry(Union{Real,AbstractVector{<:Real}}, "-", "Radial location of the deposition profile"; default=0.0)
    current_efficiency = Entry(Union{Real,AbstractVector{<:Real}}, "A/W", "Current drive efficiency"; default=0.3)
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
        gaussian_source_to_dd(
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
            2;
            electrons_particles=beam_particles,
            momentum_tor=momentum_source,
            j_parallel=j_parallel
        )
    end
    return actor
end

#= == =#
#  EC  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorECsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width = Entry(Union{Real,AbstractVector{<:Real}}, "-", "Width of the deposition profile"; default=0.1)
    rho_0 = Entry(Union{Real,AbstractVector{<:Real}}, "-", "Radial location of the deposition profile"; default=0.0)
    current_efficiency = Entry(Union{Real,AbstractVector{<:Real}}, "A/W", "Current drive efficiency"; default=0.2)
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

        isource = resize!(cs.source, "identifier.name" => ecl.name)
        gaussian_source_to_dd(
            isource,
            ecl.name,
            3,
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

#= == =#
#  IC  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorICsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width = Entry(Union{Real,AbstractVector{<:Real}}, "-", "Width of the deposition profile"; default=0.1)
    rho_0 = Entry(Union{Real,AbstractVector{<:Real}}, "-", "Radial location of the deposition profile"; default=0.0)
    current_efficiency = Entry(Union{Real,AbstractVector{<:Real}}, "A/W", "Current drive efficiency"; default=0.125)
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
        gaussian_source_to_dd(
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

#= == =#
#  LH  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorLHsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width = Entry(Union{Real,AbstractVector{<:Real}}, "-", "Width of the deposition profile"; default=0.15)
    rho_0 = Entry(Union{Real,AbstractVector{<:Real}}, "-", "Radial location of the deposition profile"; default=0.6)
    current_efficiency = Entry(Union{Real,AbstractVector{<:Real}}, "A/W", "Current drive efficiency"; default=0.4)
end

mutable struct ActorLHsimple <: HCDAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorLHsimple
    width::AbstractVector{<:Real}
    rho_0::AbstractVector{<:Real}
    current_efficiency::AbstractVector{<:Real}
end

"""
    ActorLHsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the Lower-hybrid electron energy deposition and current drive as a gaussian.

!!! note 
    Stores data in `dd.lh_antennas, dd.core_sources`
"""
function ActorLHsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorLHsimple(kw...)
    actor = ActorLHsimple(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorLHsimple(dd::IMAS.dd, par::FUSEparameters__ActorLHsimple; kw...)
    logging_actor_init(ActorLHsimple)
    par = par(kw...)
    n_antennas = length(dd.lh_antennas.antenna)
    _, width, rho_0, current_efficiency = same_length_vectors(1:n_antennas, par.width, par.rho_0, par.current_efficiency)
    return ActorLHsimple(dd, par, width, rho_0, current_efficiency)
end

function _step(actor::ActorLHsimple)
    for (idx, lha) in enumerate(actor.dd.lh_antennas.antenna)
        eqt = actor.dd.equilibrium.time_slice[]
        cp1d = actor.dd.core_profiles.profiles_1d[]
        cs = actor.dd.core_sources

        power_launched = @ddtime(lha.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
        area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

        ion_electron_fraction = 0.0

        ne_vol = integrate(volume_cp, cp1d.electrons.density) / volume_cp[end]
        j_parallel = actor.current_efficiency[idx] / eqt.boundary.geometric_axis.r / (ne_vol / 1e19) * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        isource = resize!(cs.source, "identifier.name" => lha.name)
        gaussian_source_to_dd(
            isource,
            lha.name,
            4,
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

#= ========= =#
#  functions  #
#= ========= =#
function gaussian_source_to_dd(
    isource,
    name,
    index,
    rho_cp,
    volume_cp,
    area_cp,
    power_launched,
    ion_electron_fraction,
    rho_0,
    width,
    gauss_order;
    electrons_particles=missing,
    momentum_tor=missing,
    j_parallel=missing
)
    gaussian = sgaussian(rho_cp, rho_0, width, gauss_order)
    gaussian_vol = gaussian / integrate(volume_cp, gaussian)
    gaussian_area = gaussian / integrate(area_cp, gaussian)

    electrons_energy = power_launched .* gaussian_vol .* (1 .- ion_electron_fraction)
    if sum(electrons_energy) == 0.0
        electrons_energy = missing
    end

    total_ion_energy = power_launched .* gaussian_vol .* ion_electron_fraction
    if sum(total_ion_energy) == 0.0
        total_ion_energy = missing
    end

    if electrons_particles !== missing
        electrons_particles = gaussian_vol .* electrons_particles
    end
    if momentum_tor !== missing
        momentum_tor = gaussian_area .* momentum_tor
    end
    if j_parallel !== missing
        j_parallel = gaussian_area .* j_parallel
    end
    return IMAS.new_source(isource, index, name, rho_cp, volume_cp; electrons_energy, total_ion_energy, electrons_particles, j_parallel, momentum_tor)
end

function sgaussian(rho::Union{LinRange,Vector}, rho_0::Real, width::Real, order::Real=1.0)
    return exp.(-((rho .- rho_0) .^ 2 / 2width^2) .^ order)
end
