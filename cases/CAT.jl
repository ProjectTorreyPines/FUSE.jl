function Parameters(::Type{Val{:CAT}})
    par = Parameters()
    par.general.casename = "CAT"
    par.general.init_from = :ods

    par.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "CAT_eq_ods.json")

    par.build.blanket = 1.0
    par.build.shield = 0.5
    par.build.vessel = 0.125
    par.build.symmetric = true

    par.pf_active.n_oh_coils = 6
    par.pf_active.n_pf_coils_inside = 0
    par.pf_active.n_pf_coils_outside = 6
    par.pf_active.technology = coil_technology(:ITER, :PF)

    par.tf.shape = :triple_arc
    par.tf.n_coils = 16
    par.tf.technology = coil_technology(:ITER, :TF)

    par.oh.technology = coil_technology(:ITER, :OH)
    par.oh.flattop_duration = 1000

    par.core_profiles.ne_ped = 7E19
    par.core_profiles.n_peaking = 1.5
    par.core_profiles.T_shaping = 1.8
    par.core_profiles.w_ped = 0.08
    par.core_profiles.zeff = 2.5
    par.core_profiles.rot_core = 0.0
    par.core_profiles.bulk = :DT
    par.core_profiles.impurity = :Ne

    par.nbi.beam_power = 20E6
    par.nbi.beam_energy = 200e3
    par.nbi.beam_mass = 2
    par.nbi.toroidal_angle = 0.0

    return set_new_base!(par)
end