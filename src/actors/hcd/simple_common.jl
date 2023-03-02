import NumericalIntegration: integrate

#= ========= =#
#  functions  #
#= ========= =#
function gaussian_source(
    isource::Integer,
    name::AbstractString,
    index::Integer,
    rho_cp::Real,
    volume_cp::Real,
    area_cp::Real,
    power_launched::Real,
    ion_electron_fraction::Real,
    rho_0::Real,
    width::Real,
    gauss_order::Float64;
    electrons_particles=missing,
    momentum_tor=missing,
    j_parallel=missing)

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

function sgaussian(rho::AbstractVector{<:Real}, rho_0::Real, width::Real, order::Float64=1.0)
    return exp.(-((rho .- rho_0) .^ 2 / 2width^2) .^ order)
end
