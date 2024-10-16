import FuseUtils: trapz

#= ========= =#
#  functions  #
#= ========= =#
function shaped_source(
    source::IMAS.core_sources__source,
    name::AbstractString,
    index::Integer,
    rho_cp::Vector{<:Real},
    volume_cp::Vector{<:Real},
    area_cp::Vector{<:Real},
    power_launched::Real,
    ion_electron_fraction_cp::Vector{<:Real},
    shape_function::Function;
    electrons_particles=missing,
    momentum_tor=missing,
    j_parallel=missing)

    gaussian = shape_function.(rho_cp)
    gaussian_vol = gaussian / trapz(volume_cp, gaussian)
    gaussian_area = gaussian / trapz(area_cp, gaussian)

    electrons_energy = power_launched .* gaussian_vol .* (1 .- ion_electron_fraction_cp)
    if sum(electrons_energy) == 0.0
        electrons_energy = missing
    end

    total_ion_energy = power_launched .* gaussian_vol .* ion_electron_fraction_cp
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

    return IMAS.new_source(source, index, name, rho_cp, volume_cp, area_cp; electrons_energy, total_ion_energy, electrons_particles, j_parallel, momentum_tor)
end
