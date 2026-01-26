#= ============= =#
#  Simple PELLET  #
#= ============= =#
@actor_parameters_struct _ActorSimplePLactuator{T} begin
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.5, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.25, check=x -> @assert x >= 0.0 "must be: width > 0.0")
end

@actor_parameters_struct ActorSimplePL{T} begin
    actuator::ParametersVector{_FUSEparameters__ActorSimplePLactuator{T}} = ParametersVector{_FUSEparameters__ActorSimplePLactuator{T}}()
end

mutable struct ActorSimplePL{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSimplePL{P}}
    function ActorSimplePL(dd::IMAS.dd{D}, par::FUSEparameters__ActorSimplePL{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSimplePL)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSimplePL(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates pellet fueling particle deposition using a simplified deposition model.

This actor models pellet particle injection through:
- Calculation of pellet particle content based on geometry (spherical, cylindrical, or rectangular)
- Material density lookup for different species (DT, D, T, C, Ne)
- Frequency-based particle source calculation
- Beta-function radial deposition profile with configurable location and width
- Pure particle source (no heating or current drive)

The pellet volume and particle content are calculated based on:
- Spherical: V = (4/3)πr³
- Cylindrical: V = πr²h  
- Rectangular: V = l×w×h

Particle rate = pellet_volume × material_density × injection_frequency

The deposition profile uses a beta function with parameters tuned for pellet-like deposition.

!!! note

    Reads data in `dd.pellets.launcher`, `dd.pulse_schedule` and stores results in `dd.core_sources`
"""
function ActorSimplePL(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSimplePL(dd, act.ActorSimplePL; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorSimplePL)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    cs = dd.core_sources

    rho_cp = cp1d.grid.rho_tor_norm
    volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
    area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

    for (k,(ps, pll)) in enumerate(zip(dd.pulse_schedule.pellet.launcher, dd.pellets.launcher))
        frequency = @ddtime(ps.frequency.reference)
        rho_0 = par.actuator[k].rho_0
        width = par.actuator[k].width

        shape = IMAS.index_2_name(pll.shape)[pll.shape.type.index]
        size = pll.shape.size

        density = pellet_density(pll.species[1].label)
        if shape == :spherical
            volume = (4.0 / 3.0) * π * size[1]^3
        elseif shape == :cylindrical
            volume = π * size[1]^2 * size[2]
        elseif shape == :rectangular
            volume = size[1] * size[2] * size[3]
        end
        electrons_particles = volume * density * frequency

        ion_electron_fraction_cp = zeros(length(rho_cp))

        source = resize!(cs.source, :pellet, "identifier.name" => pll.name; wipe=false)

        mode = -0.9
        α = 2.0
        β = ((α - 1) / mode) - (α - 2)
        ρ_peak = (α - 1) / (α + β - 2) / 2 + 0.5
        beta_width = width * 4

        shaped_source!(
            source,
            pll.name,
            source.identifier.index,
            rho_cp,
            volume_cp,
            area_cp,
            0.0,
            ion_electron_fraction_cp,
            ρ -> IMAS.beta(ρ + ρ_peak * beta_width, rho_0, beta_width, mode);
            electrons_particles
        )
    end

    return actor
end

# the number of particle/m^3 # the data of mass density and atomic weight is sourced from PAM model within OMFIT
function pellet_density(species::String)
    material_density = Dict("DT" => 0.257, "D" => 0.2, "T" => 0.318, "C" => 3.3, "Ne" => 1.44)
    atomic_weigth = Dict("DT" => 2.515, "D" => 2.014, "T" => 3.016, "C" => 12.011, "Ne" => 20.183)
    return material_density[species] * IMAS.mks.avog / atomic_weigth[species] * 1E6
end