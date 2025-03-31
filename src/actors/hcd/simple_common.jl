import IMASutils: trapz

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

function setup(ecb::IMAS.ec_launchers__beam, eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall, par::_FUSEparameters__ActorSimpleECactuator)
    # Estimate operating frequency and mode
    if ismissing(ecb.frequency, :data)
        resonance = IMAS.ech_resonance(eqt)
        ecb.frequency.time = [-Inf]
        ecb.frequency.data = [resonance.frequency]
        ecb.mode = resonance.mode == "X" ? -1 : 1
    end
    # Pick a reasonable launch location
    if ismissing(ecb.launching_position, :r) || ismissing(ecb.launching_position, :z)
        fw = IMAS.first_wall(wall)
        if !isempty(fw.r)
            index = argmax(fw.r .+ fw.z)
            @ddtime(ecb.launching_position.r = fw.r[index])
            @ddtime(ecb.launching_position.z = fw.z[index])
        else
            index = argmax(eqt.boundary.outline.r .+ eqt.boundary.outline.z)
            @ddtime(ecb.launching_position.r = eqt.boundary.outline.r[index])
            @ddtime(ecb.launching_position.z = eqt.boundary.outline.z[index])
        end
    end
    if ismissing(ecb.launching_position, :phi)
        @ddtime(ecb.launching_position.phi = 0.0)
    end
    # beam properties
    if ismissing(ecb.phase, :angle) || ismissing(ecb.phase, :curvature)
        @ddtime(ecb.phase.angle = 0.0)
        @ddtime(ecb.phase.curvature = [0.0, 0.0])
    end
    if ismissing(ecb.spot, :angle) || ismissing(ecb.spot, :size)
        @ddtime(ecb.spot.angle = 0.0)
        @ddtime(ecb.spot.size = [0.0172, 0.0172])
    end
    # aiming based on rho0
    if (ismissing(ecb, :steering_angle_tor) || ismissing(ecb, :steering_angle_pol)) && !ismissing(par, :rho_0)
        launch_r = @ddtime(ecb.launching_position.r)
        launch_z = @ddtime(ecb.launching_position.z)
        resonance_layer = IMAS.ech_resonance_layer(eqt, IMAS.frequency(ecb))
        _, _, RHO_interpolant = IMAS.ρ_interpolant(eqt)
        rho_resonance_layer = RHO_interpolant.(resonance_layer.r, resonance_layer.z)
        index = resonance_layer.z .> eqt.global_quantities.magnetic_axis.z
        sub_index = argmin(abs.(rho_resonance_layer[index] .- par.rho_0))
        @ddtime(ecb.steering_angle_tor = 0.0)
        @ddtime(ecb.steering_angle_pol = atan(resonance_layer.r[index][sub_index] - launch_r, resonance_layer.z[index][sub_index] - launch_z) + pi / 2)
    end
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