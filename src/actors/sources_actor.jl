import NumericalIntegration: integrate, cumul_integrate

# simple core_sources initialization which runs the simple_HCD_actors

function add_source(
    source::IMAS.core_sources__source,
    index::Int,
    name::String,
    rho::Union{Vector{},LinRange},
    volume::Union{Vector{},LinRange};
    Qe = missing,
    Qi = missing,
    qe = missing,
    qi = missing,
    se = missing,
    jpar = missing,
    momentum = missing)

    source.identifier.name = name
    source.identifier.index = index
    resize!(source.profiles_1d)
    cs1d = source.profiles_1d[]
    cs1d.grid.rho_tor_norm = rho
    cs1d.grid.volume = volume

    if qe !== missing
        cs1d.electrons.energy = IMAS.interp(LinRange(0, 1, length(qe)), qe)(cs1d.grid.rho_tor_norm)
    end
    if Qe !== missing
        cs1d.electrons.power_inside = IMAS.interp(LinRange(0, 1, length(Qe)), Qe)(cs1d.grid.rho_tor_norm)
    end
    if se !== missing
        cs1d.electrons.particles = IMAS.interp(LinRange(0, 1, length(se)), se)(cs1d.grid.rho_tor_norm)
    end
    if qi !== missing
        cs1d.total_ion_energy = IMAS.interp(LinRange(0, 1, length(qi)), qi)(cs1d.grid.rho_tor_norm)
    end
    if Qi !== missing
        cs1d.total_ion_power_inside = IMAS.interp(LinRange(0, 1, length(Qi)), Qi)(cs1d.grid.rho_tor_norm)
    end
    if momentum !== missing
        cs1d.momentum_tor = IMAS.interp(LinRange(0, 1, length(momentum)), momentum)(cs1d.grid.rho_tor_norm)
    end
    if jpar !== missing
        cs1d.jpar = IMAS.interp(LinRange(0, 1, length(jpar)), jpar)(cs1d.grid.rho_tor_norm)
    end
end

function sgaussian(rho::Union{LinRange,Vector}, rho_0::Real, width::Real, order::Real = 1.0)
    return exp.(-((rho .- rho_0) .^ 2 / 2width^2) .^ order)
end

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
    # Calculates the nbi power deposition fraction profile from sivukhin
    x = beam_energy ./ W_crit
    y = x .* rho
    f = integrate(y, vec(1.0 / (1.0 .+ y .^ 1.5)))
    return f ./ x
end

function spitzer_conductivity(ne, Te, Zeff)
    # Calculates the spitzer conductivity
    lnLambda_e = 31.3 - log(sqrt(ne) / (Te))
    NZ_term = 0.58 + 0.74 / (0.76 + Zeff)
    return 1.9012e4 * Te .^ 1.5 / (Zeff * NZ_term * lnLambda_e)
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


# make the NBI/EC/yada yda init and actor

# init core_sources will populate dd.nbi. dd.eclauchner, runs the simple_actors adn populate dd.core_sources.

# simple_nbiactor, simple_ecactor ....

# simple_nbiactor will setup from nbi.dd and run actor

# TGYRO-like solver 


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

function init_nbis(dd::IMAS.dd; beam_energy::Real, beam_mass::Real, beam_power::Real, toroidal_angle::Real)
    return init_nbi(dd::IMAS.dd; beam_energy=[beam_energy], beam_mass=[beam_mass], beam_power=[beam_power], toroidal_angle=[toroidal_angle])
end
# beam_energy=200e3, beam_mass = 2., nbi_width = 0.18, nbi_current_efficiency = 0.3, nbi_rho_0 = 0.0, injection_angle=19.6

function init_nbi(dd::IMAS.dd; beam_energy::Vector, beam_mass::Vector, beam_power::Vector, toroidal_angle::Vector)
    for idx in 1:length(beam_power)
        resize!(dd.nbi.unit,1)
        nbi_u = dd.nbi.unit[idx]
        @ddtime (nbi_u.energy.data = beam_energy[idx])
        @ddtime (nbi_u.power_launched.data = beam_power[idx])
        nbi_u.species.a = beam_mass[idx]
        # 1 beamlet
        resize!(nbi_u.beamlets_group, 1)
        nbi_u.beamlets_group[1].angle = toroidal_angle[idx] / 360 * 2pi
    end
end


function init_simple_core_sources(dd::IMAS.dd;
    power_nbi = missing, power_ec = missing, power_ic = missing,
    power_lh = missing, power_ohm = missing)

    if power_nbi !== missing
        init_nbi(dd, 200e3, 2., power_nbi, 19.6)
        nbiactor = simpleNBIactor(dd)
        FUSE.step(nbiactor)
    end
end

function step(actor::simpleNBIactor; verbose=false)
    for (idx, nbi_u) in enumerate(actor.dd.nbi.unit)
        eqt = actor.dd.equilibrium.time_slice[]
        cp1d = actor.dd.core_profiles.profiles_1d[]
        cs = actor.dd.core_sources

        beam_energy = @ddtime (nbi_u.energy.data)
        beam_mass = nbi_u.species.a
        power_launched = @ddtime(nbi_u.power_launched.data)

        # Construct nbi source
        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume)[rho_cp]
        area_cp = IMAS.interp(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area)[rho_cp]

        nbi_gaussian = sgaussian(rho_cp, actor.rho_0[idx], actor.width[idx], 2)

        qi = power_launched * sivukhin_fraction(beam_energy, cp1d) .* nbi_gaussian
        qi /= integrate(volume_cp, qi)

        qe = power_launched * (1 .- sivukhin_fraction(beam_energy, cp1d)) .* nbi_gaussian
        qe /= integrate(volume_cp, qe)

        beam_particles = power_launched / (beam_energy * constants.e)
        se = beam_particles .* nbi_gaussian
        se /= integrate(volume_cp, se)

        momentum =  sin(nbi_u.beamlets_group[1].angle) * beam_particles * sqrt(2 * beam_energy * constants.e / beam_mass / constants.m_u) * beam_mass * constants.m_u .* nbi_gaussian
        momentum /= integrate(area_cp, momentum)

        isource = resize!(cs.source, "identifier.name" => "beam_$idx")
        @show typeof(qe),typeof(volume_cp)
        add_source(isource, 8, "beam_$idx", rho_cp, volume_cp; qe=qe , Qe=cumul_integrate(volume_cp, qe), qi=qi , Qi=cumul_integrate(volume_cp, qi) , se=se , momentum=momentum)
    end
end