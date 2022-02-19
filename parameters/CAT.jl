struct CATparameters end
parametersDispatcher[:CAT] = CATparameters()

function Parameters(::CATparameters)
    params = Parameters()
    params.general.init_from = :ods

    params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "CAT_eq_ods.json")

    params.build.is_nuclear_facility = false

    params.pf_active.n_oh_coils = 6
    params.pf_active.n_pf_coils_inside = 0
    params.pf_active.n_pf_coils_outside = 6

    params.core_profiles.ne_ped = 7E19
    params.core_profiles.n_peaking = 1.5
    params.core_profiles.T_shaping = 1.8
    params.core_profiles.w_ped = 0.08
    params.core_profiles.zeff = 2.5
    params.core_profiles.rot_core = 0.0
    params.core_profiles.bulk = :DT
    params.core_profiles.impurity = :Ne

    params.nbi.beam_power = 20E6
    params.nbi.beam_energy = 200e3
    params.nbi.beam_mass = 2
    params.nbi.toroidal_angle = 0.0

    return parameters
end