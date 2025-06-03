#= ========= =#
#  functions  #
#= ========= =#
function shaped_source!(
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
        momentum_tor = gaussian_vol .* momentum_tor
    end

    if j_parallel !== missing
        j_parallel = gaussian_area .* j_parallel
    end

    return IMAS.new_source(source, index, name, rho_cp, volume_cp, area_cp; electrons_energy, total_ion_energy, electrons_particles, j_parallel, momentum_tor)
end

"""
    populate_wave1d_from_source1d!(wv1d::IMAS.waves__coherent_wave___profiles_1d, cs1d::IMAS.core_sources__source___profiles_1d)

Fills out waves.coherent_wave[:].profiles_1d[] based on the data found in core_sources.source[:].profiles_1d[]
"""
function populate_wave1d_from_source1d!(wv1d::IMAS.waves__coherent_wave___profiles_1d, cs1d::IMAS.core_sources__source___profiles_1d)
    wv1d.grid.rho_tor_norm = cs1d.grid.rho_tor_norm

    # currents
    if !ismissing(cs1d, :j_parallel)
        wv1d.current_parallel_density = cs1d.j_parallel
        wv1d.current_tor_inside = cs1d.current_parallel_inside
    end

    # electrons
    energy = zero(wv1d.grid.rho_tor_norm)
    power_inside = zero(wv1d.grid.rho_tor_norm)
    if !ismissing(cs1d.electrons, :energy)
        energy .+= cs1d.electrons.energy
        wv1d.electrons.power_density_thermal = cs1d.electrons.energy
        power_inside .+= cs1d.electrons.power_inside
        wv1d.electrons.power_inside_thermal = cs1d.electrons.power_inside
    end

    # ions
    resize!(wv1d.ion, length(cs1d.ion))
    for (wv1dion, cs1dion) in zip(wv1d.ion, cs1d.ion)
        wv1dion.label = cs1dion.label
        resize!(wv1dion.element, 1)
        wv1dion.element[1].a = cs1dion.element[1].a
        wv1dion.element[1].z_n = cs1dion.element[1].z_n
        if !ismissing(cs1dion, :energy)
            energy .+= cs1dion.energy
            wv1dion.power_density_thermal = cs1dion.energy
            power_inside .+= cs1dion.power_inside
            wv1dion.power_inside_thermal = cs1dion.power_inside
        end
    end

    # totals
    wv1d.power_density = energy
    wv1d.power_inside = power_inside

    return wv1d
end