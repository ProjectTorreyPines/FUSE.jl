function Parameters(::Type{Val{:D3D}})
    par = Parameters()
    par.general.casename = "D3D"
    par.general.init_from = :ods

    par.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "D3D_eq_ods.json")

    par.build.is_nuclear_facility = false
    par.build.symmetric = true

    par.pf_active.n_oh_coils = 10
    par.pf_active.n_pf_coils_inside = 8
    par.pf_active.n_pf_coils_outside = 0

    par.tf.shape = 3
    par.tf.n_coils = 18

    par.core_profiles.ne_ped = 5E19
    par.core_profiles.n_peaking = 1.5
    par.core_profiles.T_shaping = 1.8
    par.core_profiles.w_ped = 0.08
    par.core_profiles.zeff = 2.0
    par.core_profiles.rot_core = 5E3
    par.core_profiles.bulk = :D
    par.core_profiles.impurity = :C

    par.nbi.beam_power = 5E6
    par.nbi.beam_energy = 80e3
    par.nbi.beam_mass = 2
    par.nbi.toroidal_angle = 20.0 / 180 * pi

    return par
end