"""
    case_parameters(:KDEMO)
"""
function case_parameters(::Type{Val{:KDEMO}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ec=1, n_ic=1)
    act = ParametersActors()

    ini.general.casename = "K-DEMO"
    ini.general.init_from = :scalars

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 1.180137712438599,
        :OH => 0.35554973675585266,
        :hfs_TF => 0.3074802716013343,
        :hfs_gap_TF_vacuum_vessel => 0.13784092291547867,
        :hfs_vacuum_vessel_wall_outer => 0.02,
        :hfs_vacuum_vessel => 0.13230115364434794,
        :hfs_vacuum_vessel_wall_inner => 0.02,
        :hfs_gap_high_temp_shield_vacuum_vessel => 0.06892046145773945,
        :hfs_high_temp_shield => 0.2067613843732179,
        :hfs_blanket => 0.3316196291966129,
        :hfs_first_wall => 0.02,
        :plasma => 3.1014207655982666,
        :lfs_first_wall => 0.02,
        :lfs_blanket => 0.7654012311802618,
        :lfs_high_temp_shield => 0.2067613843732179,
        :lfs_gap_high_temp_shield_vacuum_vessel => 0.3446023072886959,
        :lfs_vacuum_vessel_wall_inner => 0.02,
        :lfs_vacuum_vessel => 0.13230115364434794,
        :lfs_vacuum_vessel_wall_outer => 0.02,
        :lfs_gap_TF_vacuum_vessel => 0.13784092291547867,
        :lfs_TF => 0.3074802716013334,
        :gap_cryostat => 1.4221989470234107,
        :cryostat => 0.21209509270084714
    )
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = false
    ini.build.divertors = :lower
    ini.build.n_first_wall_conformal_layers = 1

    ini.equilibrium.B0 = 7.5
    ini.equilibrium.R0 = 6.8
    ini.equilibrium.ϵ = 0.3088235294117647
    ini.equilibrium.κ = 1.85
    ini.equilibrium.δ = 0.485
    ini.equilibrium.pressure_core = 964500.0
    ini.equilibrium.ip = 1.2e7
    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.675
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne
    ini.core_profiles.helium_fraction = 0.01

    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 6
    ini.pf_active.technology = :nb3sn

    ini.tf.shape = :double_ellipse
    ini.tf.n_coils = 18
    # Table 2, NF 55 (2015) 053027 - KDEMO TF made of high-Jc Nb3Sn, all other coils from ITER-type Nb3Sn
    ini.tf.technology = :nb3sn_kdemo

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    ini.oh.n_coils = 6
    ini.oh.technology = :nb3sn

    ini.ec_launcher[1].power_launched = 5.0e7

    ini.ic_antenna[1].power_launched = 5.0e7

    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 1.1
    ini.requirements.power_electric_net = 400e6 # as example


    act.ActorTGLF.user_specified_model = "sat1_em_iter"
    act.ActorStabilityLimits.raise_on_breach = false

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
