import NumericalIntegration: integrate

#= ====================== =#
#     simple H&CD actors   #
#= ====================== =#
mutable struct simpleNBIactor <: AbstractActor
    dd :: IMAS.dd
    width :: Vector
    rho_0 :: Vector
    current_efficiency :: Vector
end

function simpleNBIactor(dd::IMAS.dd, width::Real=0.3, rho_0::Real = 0.0, current_efficiency::Real = 0.3)
    nbeam = ones(length(dd.nbi.unit))
    return simpleNBIactor(dd, nbeam .* width, nbeam .* rho_0, nbeam .* current_efficiency)
end

function step(actor::simpleNBIactor; verbose=false)
    for (idx, nbi_u) in enumerate(actor.dd.nbi.unit)
        eqt = actor.dd.equilibrium.time_slice[]
        cp1d = actor.dd.core_profiles.profiles_1d[]
        cs = actor.dd.core_sources

        beam_energy = @ddtime (nbi_u.energy.data)
        beam_mass = nbi_u.species.a
        power_launched = @ddtime(nbi_u.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume)[rho_cp]
        area_cp = IMAS.interp(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area)[rho_cp]

        nbi_gaussian = sgaussian(rho_cp, actor.rho_0[idx], actor.width[idx], 2)
        nbi_gaussian_vol = nbi_gaussian / integrate(volume_cp, nbi_gaussian)
        nbi_gaussian_area = nbi_gaussian / integrate(area_cp, nbi_gaussian)

        ion_elec_ratio = sivukhin_fraction(beam_energy, cp1d)

        total_ion_energy = power_launched .* nbi_gaussian_vol .* ion_elec_ratio
        electrons_energy = power_launched .* nbi_gaussian_vol .* (1 .- ion_elec_ratio)

        beam_particles = power_launched / (beam_energy * constants.e)
        electrons_particles = beam_particles .* nbi_gaussian_vol

        ne_vol = integrate(volume_cp,cp1d.electrons.density) / volume_cp[end]
        j_parallel = nbi_gaussian_area .* actor.current_efficiency / eqt.boundary.geometric_axis.r / (ne_vol/1e19) * power_launched
        momentum_source= sin(nbi_u.beamlets_group[1].angle) * beam_particles * sqrt(2 * beam_energy * constants.e / beam_mass / constants.m_u) * beam_mass * constants.m_u 
        momentum_tor =  nbi_gaussian_area .* momentum_source

        isource = resize!(cs.source, "identifier.name" => "beam_$idx")
        IMAS.new_source(isource, 8, "beam_$idx", rho_cp, volume_cp; electrons_energy, total_ion_energy, electrons_particles, j_parallel, momentum_tor)
    end
end




"""
# Widths :
minors_radius_DIIID = 0.606
if 'aspect' in zero_d_list:
    minors_radius = zero_d_list['R'] / zero_d_list['aspect']
else:
    minors_radius = zero_d_list['a']
w0_IC = 0.15 * minors_radius_DIIID / minors_radius  # tok linear size scaling of IC and EC resonance
w0_EC = 0.15 * minors_radius_DIIID / minors_radius  # tok linear size scaling of IC and EC resonance
w0_NBI = 0.3
w0_LH = 0.1
w0_OHM = 0.0

### Current drive
# SI & KESSEL norm [W/A/m^2] / 1e19
# A DIIID case has been used to calculate gfw_dict current drive efficiencies, for calibration see OMFIT project: /fusion/projects/omfit-results/slendebroekt/projects/CHEF_PROFILE_COOKING_MARCH2.zip
# Pythontask: OMFIT['CHEF_1']['currentdrive_efficiency_calibration']
gfw_dict = {
    "NBI": 0.18,
    "OHM": 0.0,  # Ohmic currently set to 0 since ohmic calculated from total
    "EC": 0.2,
    "IC": 0.125,
    "LH": 0.43,
}  # LH current drive effciency from friedberg 2015 paper
"""
function sgaussian(rho::Union{LinRange,Vector}, rho_0::Real, width::Real, order::Real = 1.0)
    return exp.(-((rho .- rho_0) .^ 2 / 2width^2) .^ order)
end

"""
    Calculates the NBI power deposition fraction profile from sivukhin
"""
function sivukhin_fraction(beam_energy, cp1d)
    Te = cp1d.electrons.temperature
    ne = cp1d.electrons.density
    rho = cp1d.grid.rho_tor_norm

    W_crit = zeros(length(rho))
    for ion in cp1d.ion
        ni = ion.density
        Zi = ion.element[1].z_n
        mi = ion.element[1].a    
        W_crit .+= Te .* 14.8 .* (mi^1.5 .* ni ./ ne .* Zi .* Zi ./ mi) .^ 0.666  # [P]
    end

    x = beam_energy ./ W_crit
    y = x .* rho
    f = integrate(y, vec(1.0 / (1.0 .+ y .^ 1.5)))
    ion_elec_fraction = f ./ x

    if any(i -> i > 1.0, ion_elec_fraction)
        error("NBI power is going more to ions than electrons. Check that beam_energy $(beam_energy) in eV")
    end

    return ion_elec_fraction
end