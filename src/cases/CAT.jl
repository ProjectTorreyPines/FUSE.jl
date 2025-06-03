"""
    case_parameters(:CAT)

GA Compact Advanced Tokamak design
"""
function case_parameters(::Type{Val{:CAT}})
    ini = ParametersInits()
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
        vessel=0.5,
        pf_inside_tf=false,
        pf_outside_tf=true,
        thin_vessel_walls=true)

    ini.build.n_first_wall_conformal_layers = 1

    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:gap_cryostat].coils_inside = 6

    ini.oh.technology = :nb3sn_iter
    ini.pf_active.technology = :nb3sn_iter
    ini.tf.technology = :nb3sn_iter

    ini.tf.shape = :circle_ellipse
    ini.tf.n_coils = 16

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.8 * 0.75
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.zeff = 2.5
    ini.core_profiles.helium_fraction = 0.01
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    ini.requirements.flattop_duration = 1000.0
    ini.requirements.tritium_breeding_ratio = 1.1

    resize!(ini.nb_unit, 1)
    ini.nb_unit[1].power_launched = 20E6
    ini.nb_unit[1].beam_energy = 200e3
    ini.nb_unit[1].beam_mass = 2.0
    ini.nb_unit[1].normalized_tangency_radius = 0.8
    ini.nb_unit[1].beam_current_fraction = [1.0,0.0,0.0]
    ini.nb_unit[1].current_direction = :co
    ini.nb_unit[1].offaxis = false

    #### ACT ####

    act.ActorPFdesign.symmetric = true

    act.ActorFluxMatcher.algorithm = :simple

    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    return ini, act
end
