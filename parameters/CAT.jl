struct CATparameters end
parametersDispatcher[:CAT] = CATparameters()

function Parameters(::CATparameters)
    par = Parameters()
    par.general.init_from = :ods

    par.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "CAT_eq_ods.json")

    par.build.is_nuclear_facility = false

    par.pf_active.n_oh_coils = 6
    par.pf_active.n_pf_coils_inside = 0
    par.pf_active.n_pf_coils_outside = 6

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

    return parameters
end