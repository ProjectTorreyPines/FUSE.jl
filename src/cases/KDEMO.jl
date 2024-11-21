"""
    case_parameters(:KDEMO)
"""
function case_parameters(::Type{Val{:KDEMO}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()

    ini.general.casename = "K-DEMO"
    ini.general.init_from = :scalars

    gaps_thickness = 0.1
    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 1.597,
        :OH => 0.554,
        :hfs_TF => 1.204,
        :hfs_gap_TF_vacuum_vessel => gaps_thickness,
        :hfs_vacuum_vessel_outer => 0.1,
        :hfs_gap_water => 0.15,
        :hfs_vacuum_vessel_inner => 0.1,
        :hfs_gap_high_temp_shield_vacuum_vessel => 0.01,
        :hfs_high_temp_shield => 0.2,
        :hfs_blanket => 0.38,
        :hfs_first_wall => 0.02,
        :plasma => 3.1,
        :lfs_first_wall => 0.02,
        :lfs_blanket => 1.2,
        :lfs_high_temp_shield => 0.2,
        :lfs_gap_high_temp_shield_vacuum_vessel => 0.01,
        :lfs_vacuum_vessel_inner => 0.1,
        :lfs_gap_water => 0.15,
        :lfs_vacuum_vessel_outer => 0.1,
        :lfs_gap_TF_vacuum_vessel => 0.789,
        :lfs_TF => 1.204,
        :gap_cryostat => 2.0,
        :cryostat => 0.1,
        :gap_world => 1.0,
    )

    ini.build.plasma_gap = 0.125
    ini.build.symmetric = false
    ini.build.divertors = :lower
    ini.build.n_first_wall_conformal_layers = 1

    ini.equilibrium.B0 = 6.55
    ini.equilibrium.R0 = 6.8
    ini.equilibrium.ϵ = 0.3088235294117647
    ini.equilibrium.κ = 2.0
    ini.equilibrium.δ = 0.59
    ini.equilibrium.ip = 1.3e7
    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.675
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 25E3
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne
    ini.core_profiles.helium_fraction = 0.01

    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:gap_cryostat].coils_inside = 6

    ini.oh.technology = :nb3sn
    ini.pf_active.technology = :nb3sn
    # Table 2, NF 55 (2015) 053027 - KDEMO TF made of high-Jc Nb3Sn, all other coils from ITER-type Nb3Sn
    ini.tf.technology = :nb3sn_kdemo

    ini.tf.shape = :double_ellipse
    ini.tf.n_coils = 18

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    resize!(ini.ec_launcher, 1)
    ini.ec_launcher[1].power_launched = 5.0e7

    resize!(ini.ic_antenna, 1)
    ini.ic_antenna[1].power_launched = 5.0e7

    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 1.1
    ini.requirements.power_electric_net = 300e6 # as example

    #### ACT ####

    act.ActorTGLF.tglfnn_model = "sat1_em_iter"

    act.ActorStabilityLimits.raise_on_breach = false

    return ini, act
end


function case_parameters(::Type{Val{:KDEMO_compact}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()

    ini.general.casename = "K-DEMO_compact"
    ini.general.init_from = :scalars

    gaps_thickness = 0.1
    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 0.8,
        :OH => 0.22,
        :hfs_TF => 0.949,
        :hfs_gap_TF_vacuum_vessel => gaps_thickness,
        :hfs_vacuum_vessel_outer => 0.1,
        :hfs_gap_water => 0.15,
        :hfs_vacuum_vessel_inner => 0.1,
        :hfs_gap_high_temp_shield_vacuum_vessel => 0.01,
        :hfs_high_temp_shield => 0.2,
        :hfs_blanket => 0.38,
        :hfs_first_wall => 0.02,
        :plasma => 1.3, 
        :lfs_first_wall => 0.02,
        :lfs_blanket => 0.7,
        :lfs_high_temp_shield => 0.2,
        :lfs_gap_high_temp_shield_vacuum_vessel => 0.01,
        :lfs_vacuum_vessel_inner => 0.1,
        :lfs_gap_water => 0.15,
        :lfs_vacuum_vessel_outer => 0.1,
        :lfs_gap_TF_vacuum_vessel => 0.789,
        :lfs_TF => 0.949,  
        :gap_cryostat => 1.0,
        :cryostat => 0.1
    )
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = false
    ini.build.divertors = :lower
    ini.build.n_first_wall_conformal_layers = 1
    ini.equilibrium.B0 = 8.35        # Updated from bt
    ini.equilibrium.R0 = 3.516       # Updated from rmajor
    ini.equilibrium.ϵ = 0.3088235294117647
    ini.equilibrium.κ = 1.8          # Updated from kappa
    ini.equilibrium.δ = 0.6          # Updated from triang
    ini.equilibrium.ip = 7.4713e6    # Updated from Ip
    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.675
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 20E3
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.helium_fraction = 0.01
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:gap_cryostat].coils_inside = 6

    ini.oh.technology = :rebco
    ini.pf_active.technology = :rebco
    ini.tf.technology = :rebco

    ini.tf.shape = :miller
    ini.tf.n_coils = 18

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    resize!(ini.ec_launcher, 1)
    ini.ec_launcher[1].power_launched = 5.0e7

    resize!(ini.ic_antenna, 1)
    ini.ic_antenna[1].power_launched = 5.0e7

    ini.requirements.flattop_duration = 3000.0
    ini.requirements.tritium_breeding_ratio = 0.8
    ini.requirements.power_electric_net = 100e6 # as example

    #### ACT ####

    act.ActorTGLF.tglfnn_model = "sat1_em_iter"

    act.ActorStabilityLimits.raise_on_breach = false

    return ini, act
end