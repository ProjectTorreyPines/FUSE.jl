"""
    case_parameters(:CAT)

GA Compact Advanced Tokamak design
"""
function case_parameters(::Type{Val{:CAT}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()

    ini.general.casename = "CAT"
    ini.general.init_from = :ods
    ini.equilibrium.boundary_from = :ods

    ini.ods.filename = joinpath(@__DIR__, "..", "sample", "CAT_eq_ods.json")

    ini.build.blanket = 1.0
    ini.build.shield = 0.5
    ini.build.vessel = 0.125
    ini.build.n_first_wall_conformal_layers = 2
    ini.material.shield = "Tungsten"
    ini.material.blanket = "FLiBe"

    ini.oh.n_coils = 6
    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 6
    ini.pf_active.technology = :ITER

    ini.tf.shape = :double_ellipse
    ini.tf.n_coils = 16
    ini.tf.technology = :ITER

    ini.oh.technology = :ITER
    ini.core_profiles.greenwald_fraction = 0.8
    ini.core_profiles.greenwald_fraction_ped = ini.core_profiles.greenwald_fraction * 0.75
    ini.core_profiles.helium_fraction = 0.01
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.zeff = 2.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    ini.requirements.flattop_duration = 1000.0

    ini.nbi.power_launched = 20E6
    ini.nbi.beam_energy = 200e3
    ini.nbi.beam_mass = 2.0
    ini.nbi.toroidal_angle = 0.0

    act.ActorPFcoilsOpt.symmetric = true

    set_new_base!(ini)
    set_new_base!(act)
    
    return ini, act
end