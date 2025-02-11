"""
    case_parameters(::Type{Val{:EXCITE}})::Tuple{ParametersAllInits,ParametersAllActors}

GA EXCITE design
"""
function case_parameters(::Type{Val{:EXCITE}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()

    #### INI ####

    ini.general.casename = "EXCITE"
    ini.general.init_from = :scalars

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 0.253,
        :OH => 0.154,
        :gap_TF_OH => 0.01,
        :hfs_TF => 0.212,
        :hfs_gap_low_temp_shield_TF => 0.01,
        :hfs_low_temp_shield => 0.101,
        :hfs_gap_vacuum_vessel_low_temp_shield => 0.01,
        :hfs_vacuum_vessel_outer => 0.05,
        :hfs_gap_water => 0.05,
        :hfs_vacuum_vessel_inner => 0.05,
        :hfs_first_wall => 0.005,
        :plasma => 1.10,
        :lfs_first_wall => 0.005,
        :lfs_vacuum_vessel_inner => 0.05,
        :lfs_gap_water => 0.05,
        :lfs_vacuum_vessel_outer => 0.05,
        :lfs_gap_vacuum_vessel_low_temp_shield => 0.01,
        :lfs_low_temp_shield => 0.10,
        :lfs_gap_low_temp_shield_TF => 0.4,
        :lfs_TF => 0.212,
        :gap_cryostat => 0.5,
        :cryostat => 0.05)

    ini.build.plasma_gap = 0.05
    ini.build.divertors = :double
    ini.build.n_first_wall_conformal_layers = 2

    ini.equilibrium.B0 = 6.0
    ini.equilibrium.R0 = 1.5
    ini.equilibrium.ϵ = 0.5 / 1.5
    ini.equilibrium.κ = 0.85
    ini.equilibrium.δ = 0.7
    ini.equilibrium.ζ = 0.05
    ini.equilibrium.𝚶 = 0.0
    ini.equilibrium.pressure_core = 1.7e6
    ini.equilibrium.ip = 6.0e6
    ini.equilibrium.xpoints = :double
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.25
    ini.core_profiles.ne_shaping = 1.5
    ini.core_profiles.Te_shaping = 2.5
    ini.core_profiles.Te_sep = 100.0
    ini.core_profiles.Ti_Te_ratio = 0.825
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.helium_fraction = 0.0
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :Kr
    ini.core_profiles.rot_core = 0.0

    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:gap_cryostat].coils_inside = 6

    ini.tf.technology = :rebco
    ini.oh.technology = :rebco
    ini.pf_active.technology = :nb3sn

    ini.tf.n_coils = 16
    ini.tf.shape = :double_ellipse

    resize!(ini.ec_launcher, 2)
    resize!(act.ActorSimpleEC.actuator, 2)
    act.ActorSimpleEC.actuator[1].rho_0 = 0.5
    act.ActorSimpleEC.actuator[1].width = 0.1
    ini.ec_launcher[1].power_launched = 25.0e6
    ini.ec_launcher[1].efficiency_conversion = 0.45
    ini.ec_launcher[1].efficiency_transmission = 0.8
    act.ActorSimpleEC.actuator[2].rho_0 = 0.7
    act.ActorSimpleEC.actuator[2].width = 0.1
    ini.ec_launcher[2].power_launched = 25.0e6
    ini.ec_launcher[2].efficiency_conversion = 0.45
    ini.ec_launcher[2].efficiency_transmission = 0.8

    ini.requirements.power_electric_net = 0.0e0
    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 0.0
    ini.requirements.coil_j_margin = 0.1
    ini.requirements.coil_stress_margin = 0.1

    ini.time.simulation_start = 0.0

    #### ACT ####

    act.ActorPlasmaLimits.models = [:q95_gt_2, :κ_controllability]

    act.ActorFluxSwing.operate_oh_at_j_crit = true

    act.ActorEquilibrium.model = :TEQUILA
    act.ActorTEQUILA.relax = 0.25

    act.ActorCoreTransport.model = :none # No flux matching possible

    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    return ini, act
end
