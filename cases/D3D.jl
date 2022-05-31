"""
    case_parameters(:D3D)

DIII-D
"""
function case_parameters(::Type{Val{:D3D}})::Tuple{ParametersInit, ParametersActor}
    ini = ParametersInit()
    act = ParametersActor()

    ini.general.casename = "D3D"
    ini.general.init_from = :ods

    ini.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "D3D_eq_ods.json")

    ini.build.blanket = 0.0
    ini.build.shield = 0.0
    ini.build.vessel = 0.0
    ini.build.n_first_wall_conformal_layers = 2

    ini.pf_active.n_oh_coils = 10
    ini.pf_active.n_pf_coils_inside = 8
    ini.pf_active.n_pf_coils_outside = 0
    ini.pf_active.technology = coil_technology(:copper)

    ini.tf.shape = :triple_arc
    ini.tf.n_coils = 24
    ini.tf.technology = coil_technology(:copper)

    ini.oh.technology = coil_technology(:copper)
    ini.oh.flattop_duration = 5

    ini.core_profiles.ne_ped = 5E19
    ini.core_profiles.greenwald_fraction = 0.7
    ini.core_profiles.helium_fraction = 0.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.w_ped = 0.08
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 5E3
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C

    ini.nbi.power_launched = 5E6
    ini.nbi.beam_energy = 80e3
    ini.nbi.beam_mass = 2
    ini.nbi.toroidal_angle = 20.0 / 180 * pi

    act.ActorPFcoilsOpt.symmetric = true

    return set_new_base!(ini), set_new_base!(act)
end