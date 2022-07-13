import NumericalIntegration: integrate

#= === =#
#  NBI  #
#= === =#
Base.@kwdef mutable struct ActorNBIsimple <: HCDAbstractActor
    dd::IMAS.dd
    width::Vector{Real}
    rho_0::Vector{Real}
    current_efficiency::Vector{Real}
end

function ParametersActor(::Type{Val{:ActorNBIsimple}})
    par = ParametersActor(nothing)
    par.width = Entry(Real, "", "Width of the deposition profile"; default=0.3)
    par.rho_0 = Entry(Real, "", "Radial location of the deposition profile"; default=0.0)
    par.current_efficiency = Entry(Real, "A/W", "Current drive efficiency"; default=0.3)
    return par
end

"""
    ActorNBIsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor estimates the NBI ion/electron energy deposition, particle source, rotation and current drive source with a super-gaussian.

!!! note 
    Stores data in `dd.nbi, dd.core_sources`
"""
function ActorNBIsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorNBIsimple(kw...)
    actor = ActorNBIsimple(dd; par.width, par.rho_0, par.current_efficiency)
    step(actor)
    finalize(actor)
    return actor
end

function ActorNBIsimple(dd::IMAS.dd; width::Real=0.3, rho_0::Real=0.0, current_efficiency::Real=0.3)
    nbeam = ones(length(dd.nbi.unit))
    return ActorNBIsimple(dd, nbeam .* width, nbeam .* rho_0, nbeam .* current_efficiency)
end

function step(actor::ActorNBIsimple)
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
Base.@kwdef mutable struct ActorECsimple <: HCDAbstractActor
    dd::IMAS.dd
    width::Vector{Real}
    rho_0::Vector{Real}
    current_efficiency::Vector{Real}
end

function ParametersActor(::Type{Val{:ActorECsimple}})
    par = ParametersActor(nothing)
    par.width = Entry(Real, "", "Width of the deposition profile"; default=0.1)
    par.rho_0 = Entry(Real, "", "Radial location of the deposition profile"; default=0.0)
    par.current_efficiency = Entry(Real, "A/W", "Current drive efficiency"; default=0.2)
    return par
end

"""
    ActorECsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor estimates the EC electron energy deposition and current drive as a gaussian.

!!! note 
    Stores data in `dd.ec_launchers, dd.core_sources`
"""
function ActorECsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorECsimple(kw...)
    actor = ActorECsimple(dd; par.width, par.rho_0, par.current_efficiency)
    step(actor)
    finalize(actor)
    return actor
end

function ActorECsimple(dd::IMAS.dd; width::Real=0.1, rho_0::Real=0.0, current_efficiency::Real=0.2)
    n_launchers = ones(length(dd.ec_launchers.launcher))
    return ActorECsimple(dd, n_launchers .* width, n_launchers .* rho_0, n_launchers .* current_efficiency)
end

function step(actor::ActorECsimple)
    for (idx, ecl) in enumerate(actor.dd.ec_launchers.launcher)
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
Base.@kwdef mutable struct ActorICsimple <: HCDAbstractActor
    dd::IMAS.dd
    width::Vector{Real}
    rho_0::Vector{Real}
    current_efficiency::Vector{Real}
end

function ParametersActor(::Type{Val{:ActorICsimple}})
    par = ParametersActor(nothing)
    par.width = Entry(Real, "", "Width of the deposition profile"; default=0.1)
    par.rho_0 = Entry(Real, "", "Radial location of the deposition profile"; default=0.0)
    par.current_efficiency = Entry(Real, "A/W", "Current drive efficiency"; default=0.125)
    return par
end

"""
    ActorICsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor estimates the ion-cyclotron electron/ion energy deposition and current drive as a gaussian.

!!! note 
    Stores data in `dd.ic_antennas, dd.core_sources`
"""
function ActorICsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorICsimple(kw...)
    actor = ActorICsimple(dd; par.width, par.rho_0, par.current_efficiency)
    step(actor)
    finalize(actor)
    return actor
end

function ActorICsimple(dd::IMAS.dd; width::Real=0.1, rho_0::Real=0.0, current_efficiency::Real=0.125)
    n_antennas = ones(length(dd.ic_antennas.antenna))
    return ActorICsimple(dd, n_antennas .* width, n_antennas .* rho_0, n_antennas .* current_efficiency)
end

function step(actor::ActorICsimple)
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
Base.@kwdef mutable struct ActorLHsimple <: HCDAbstractActor
    dd::IMAS.dd
    width::Vector{Real}
    rho_0::Vector{Real}
    current_efficiency::Vector{Real}
end

function ParametersActor(::Type{Val{:ActorLHsimple}})
    par = ParametersActor(nothing)
    par.width = Entry(Real, "", "Width of the deposition profile"; default=0.15)
    par.rho_0 = Entry(Real, "", "Radial location of the deposition profile"; default=0.6)
    par.current_efficiency = Entry(Real, "A/W", "Current drive efficiency"; default=0.4)
    return par
end

"""
    ActorLHsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor estimates the Lower-hybrid electron energy deposition and current drive as a gaussian.

!!! note 
    Stores data in `dd.lh_antennas, dd.core_sources`
"""
function ActorLHsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorLHsimple(kw...)
    actor = ActorLHsimple(dd; par.width, par.rho_0, par.current_efficiency)
    step(actor)
    finalize(actor)
    return actor
end

function ActorLHsimple(dd::IMAS.dd; width::Real=0.15, rho_0::Real=0.6, current_efficiency::Real=0.4)
    n_antennas = ones(length(dd.lh_antennas.antenna))
    return ActorICsimple(dd, n_antennas .* width, n_antennas .* rho_0, n_antennas .* current_efficiency)
end

function step(actor::ActorLHsimple)
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
