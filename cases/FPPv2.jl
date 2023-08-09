"""
    case_parameters(::Type{Val{:FPPv2}})::Tuple{ParametersAllInits,ParametersAllActors}
"""
function case_parameters(::Type{Val{:FPPv2}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()

    #### INI ####
    ini.general.casename = "FPP_v2"
    ini.general.init_from = :scalars

    ini.build.layers = OrderedCollections.OrderedDict(:gap_OH => 1.180137712438599, :OH => 0.35554973675585266, :gap_TF_OH => 0.068920461457739, :hfs_TF => 0.3074802716013343, :hfs_gap_low_temp_shield_TF => 0.06892046145773922, :hfs_low_temp_shield => 0.4241901854016936, :hfs_gap_vacuum_vessel_low_temp_shield => 0.06892046145773945, :hfs_vacuum_vessel => 0.13230115364434794, :hfs_gap_high_temp_shield_vacuum_vessel => 0.06892046145773945, :hfs_high_temp_shield => 0.2067613843732179, :hfs_blanket => 0.3316196291966129, :hfs_first_wall => 0.02, :plasma => 3.1014207655982666, :lfs_first_wall => 0.02, :lfs_blanket => 0.7654012311802618, :lfs_high_temp_shield => 0.2067613843732179, :lfs_gap_high_temp_shield_vacuum_vessel => 0.3446023072886959, :lfs_vacuum_vessel => 0.13230115364434794, :lfs_gap_vacuum_vessel_low_temp_shield => 0.068920461457739, :lfs_low_temp_shield => 0.4241901854016943, :lfs_gap_low_temp_shield_TF => 0.17230115364434795, :lfs_TF => 0.3074802716013334, :gap_cryostat => 1.4221989470234107, :cryostat => 0.21209509270084714)
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = false
    ini.build.divertors = :lower
    ini.build.n_first_wall_conformal_layers = 2

    ini.material.wall = "Tungsten"
    ini.material.blanket = "lithium-lead"
    ini.material.shield = "Steel, Stainless 316"

    ini.equilibrium.B0 = 4.713171689711136
    ini.equilibrium.R0 = 4.824432302041749
    ini.equilibrium.ϵ = 0.2857142857142857
    ini.equilibrium.κ = 0.8826
    ini.equilibrium.δ = 0.7
    ini.equilibrium.pressure_core = 1.2e6
    ini.equilibrium.ip = 8.0e6
    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.greenwald_fraction = 1.0
    ini.core_profiles.greenwald_fraction_ped = 0.7
    ini.core_profiles.T_ratio = 0.825
    ini.core_profiles.T_shaping = 2.5
    ini.core_profiles.n_shaping = 2.5
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Kr
    ini.core_profiles.helium_fraction = 0.04
    ini.core_profiles.ejima = 0.4

    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 8
    ini.pf_active.technology = :LTS

    ini.tf.n_coils = 16
    ini.tf.technology = :HTS

    ini.oh.n_coils = 6
    ini.oh.technology = :HTS

    ini.nbi.power_launched = Float64[]
    ini.nbi.beam_energy = Float64[]
    ini.nbi.efficiency_conversion = 0.4
    ini.nbi.efficiency_transmission = 0.94

    ini.ec_launchers.power_launched = 2.5e7
    ini.ec_launchers.efficiency_conversion = 0.45
    ini.ec_launchers.efficiency_transmission = 0.8

    ini.ic_antennas.power_launched = Float64[]
    ini.ic_antennas.efficiency_conversion = 0.9
    ini.ic_antennas.efficiency_transmission = 0.7
    ini.ic_antennas.efficiency_coupling = 1.0

    ini.lh_antennas.power_launched = Float64[]
    ini.lh_antennas.efficiency_conversion = 0.5
    ini.lh_antennas.efficiency_transmission = 0.72
    ini.lh_antennas.efficiency_coupling = 1.0

    ini.requirements.power_electric_net = 2.0e8
    ini.requirements.flattop_duration = 86400.0
    ini.requirements.tritium_breeding_ratio = 1.1

    #### ACT ####
    act.ActorCosting.model = :GASC
    act.ActorCosting.production_increase = 10.0
    act.ActorCosting.learning_rate = 0.616

    act.ActorECsimple.rho_0 = 0.6

    act.ActorFixedProfiles.T_shaping = 2.5
    act.ActorFixedProfiles.n_shaping = 2.5
    act.ActorFixedProfiles.T_ratio_core = 0.825

    act.ActorFluxMatcher.evolve_densities = Dict(:He4 => :match_ne_scale, :Ar => :match_ne_scale, :DT => :quasi_neutrality, :He4_fast => :constant, :electrons => :flux_match)

    act.ActorFluxSwing.operate_oh_at_j_crit = false

    act.ActorHFSsizing.j_tolerance = 0.5
    act.ActorHFSsizing.stress_tolerance = 0.01
    act.ActorHFSsizing.error_on_technology = false

    act.ActorPFcoilsOpt.optimization_scheme = :none

    act.ActorPedestal.update_core_profiles = false

    act.ActorPowerNeeds.model = :thermal_power_fraction
    act.ActorPowerNeeds.thermal_power_fraction = 0.2

    act.ActorStabilityLimits.models = [:beta_troyon_1984, :model_201, :model_401]
    act.ActorStabilityLimits.raise_on_breach = false

    act.ActorStationaryPlasma.do_plot = true
    act.ActorStationaryPlasma.convergence_error = 0.01

    act.ActorTEQUILA.relax = 0.25

    act.ActorThermalCycle.power_cycle_type = :fixed_cycle_efficiency
    act.ActorThermalCycle.fixed_cycle_efficiency = 0.4

    # finalize
    consistent_ini_act!(ini, act)
    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
