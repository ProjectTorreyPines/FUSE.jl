#= == == == =#
#  PELLET   #
#= == == == =#
Base.@kwdef mutable struct FUSEparameters__ActorPelletsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width::Entry{Union{T,AbstractVector{T}}} = Entry{Union{T,AbstractVector{T}}}("-", "Width of the deposition profile"; default=0.05)
    rho_0::Entry{Union{T,AbstractVector{T}}} = Entry{Union{T,AbstractVector{T}}}("-", "Radial location of the deposition profile"; default=0.5)
    ηcd_scale::Entry{Union{T,AbstractVector{T}}} = Entry{Union{T,AbstractVector{T}}}("-", "Scaling factor for nominal current drive efficiency"; default=1.0)
end

mutable struct ActorPelletsimple{D,P} <: HCDAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPelletsimple{P}
    function ActorPelletsimple(dd::IMAS.dd{D}, par::FUSEparameters__ActorPelletsimple{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorPelletsimple)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorPelletsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the Pellet particle deposition as a gaussian.

!!! note

    Reads data in `dd.pellet_launchers` and stores data in `dd.core_sources`
"""
function ActorPelletsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPelletsimple(dd, act.ActorPelletsimple; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorPelletsimple)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    cs = dd.core_sources

    R0 = eqt.boundary.geometric_axis.r
    rho_cp = cp1d.grid.rho_tor_norm
    volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
    area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

    #the number of particle/cm^3 # the data of mass density and atomic weight is sourced from PAM model within OMFIT
    function density_of_pellet(species::String)
        avog=6.022e23
        material_density=Dict("DT" => 0.257, "D" => 0.2, "T" => 0.318, "C" => 3.3, "Ne" => 1.44)

        atomic_weigth=Dict("DT" => 2.515, "D" => 2.014, "T" => 3.016, "C" => 12.011, "Ne" => 20.183)

        if species == "DT"
            return material_density["DT"]*avog/atomic_weigth["DT"]
        elseif species == "D"
            return material_density["D"]*avog/atomic_weigth["D"]
        elseif species == "T"
            return material_density["T"]*avog/atomic_weigth["T"]
        elseif species == "C"
            return material_density["C"]*avog/atomic_weigth["C"]
        elseif species == "Ne"
            return material_density["Ne"]*avog/atomic_weigth["Ne"]
        end
        # return is atoms / cubic centimeter 
    end
    
    _, width, rho_0 = same_length_vectors(1:length(dd.pellets.time_slice[].pellet), par.width, par.rho_0)

    for (idx, peln) in enumerate(dd.pellets.time_slice[].pellet)
        frequency=@ddtime(dd.pulse_schedule.pellets.launchers[idx].frequency.reference.data) 
        
        # Can't use @ddtime on Vector{Vector{Float64}} yet we should implement this  https://github.com/ProjectTorreyPines/IMAS.jl/issues/39
        size_matrix=dd.pulse_schedule.pellets.launchers[idx].size.reference.data
        time = dd.pulse_schedule.pellets.launchers[idx].size.reference.time
        size = size_matrix[argmin(abs.(time .- dd.global_time)),:]

        shape_array=dd.pulse_schedule.pellets.launchers[idx].shape.reference.data
        shape = shape_array[argmin(abs.(time .- dd.global_time))]

        species_array=dd.pulse_schedule.pellets.launchers[idx].species.label.reference.data
        species = species_array[argmin(abs.(time .- dd.global_time))]
        density = density_of_pellet(species)
        if shape == "spherical"
            electrons_particles=density*(4/3)pi*size[1]^3
        elseif shape == "cylinderical"
            electrons_particles=density*pi*size[1]^2*size[2]
        else shape == "rectangular"
            electrons_particles=density*size[1]*size[2]*size[3]
        end
        
        electrons_particles *= frequency
        
        ion_electron_fraction_cp = zeros(length(rho_cp))

        pellet_name = "pellet_$idx"

        source = resize!(cs.source, :pellet, "identifier.name" => pellet_name; wipe=false)
        gaussian_source(
            source,
            pellet_name,
            source.identifier.index,
            rho_cp,
            volume_cp,
            area_cp,
            0.0,
            ion_electron_fraction_cp,
            rho_0[idx],
            width[idx],
            1.0;
            electrons_particles=electrons_particles
        )
    end
    return actor
end