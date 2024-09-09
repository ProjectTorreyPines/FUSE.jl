"""
    case_parameters(:CAT)

GA Compact Advanced Tokamak design
"""
function case_parameters(::Type{Val{:CAT}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_nb=1)
    act = ParametersActors()

    ini.general.casename = "CAT"
    ini.general.init_from = :ods
    ini.equilibrium.boundary_from = :ods

    ini.ods.filename = joinpath("__FUSE__", "sample", "CAT_eq_ods.json")
    ini.time.simulation_start = 0.006

    ini.build.layers = layers_meters_from_fractions(;
        lfs_multiplier=1.0,
        wall=0.1,
        blanket=1.0,
        shield=0.5,
        vessel=0.25,
        pf_inside_tf=false,
        pf_outside_tf=true,
        thin_vessel_walls=true)
    ini.build.n_first_wall_conformal_layers = 1

    ini.oh.n_coils = 6
    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 6
    ini.pf_active.technology = :nb3sn_iter

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    ini.tf.shape = :circle_ellipse
    ini.tf.n_coils = 16
    ini.tf.technology = :nb3sn_iter

    ini.oh.technology = :nb3sn_iter
    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.8 * 0.75
    ini.core_profiles.helium_fraction = 0.01
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    ini.requirements.flattop_duration = 1000.0
    ini.requirements.tritium_breeding_ratio = 1.1

    ini.nb_unit[1].power_launched = 20E6
    ini.nb_unit[1].beam_energy = 200e3
    ini.nb_unit[1].beam_mass = 2.0
    ini.nb_unit[1].toroidal_angle = 0.0

    act.ActorPFdesign.symmetric = true

    act.ActorStabilityLimits.raise_on_breach = false

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end