struct D3Dparameters end
parametersDispatcher[:D3D] = D3Dparameters()

function Parameters(::D3Dparameters)
    params = Parameters()
    params.general.init_from = :ods

    params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "D3D_eq_ods.json")

    params.build.is_nuclear_facility = false

    params.pf_active.n_oh_coils = 20
    params.pf_active.n_pf_coils_inside = 18
    params.pf_active.n_pf_coils_outside = 0

    params.core_profiles.ne_ped = 5E19
    params.core_profiles.n_peaking = 1.5
    params.core_profiles.T_shaping = 1.8
    params.core_profiles.w_ped = 0.08
    params.core_profiles.zeff = 2.0
    params.core_profiles.rot_core = 5E3
    params.core_profiles.bulk = :D
    params.core_profiles.impurity = :C

    params.nbi.beam_power = 5E6
    params.nbi.beam_energy = 80e3
    params.nbi.beam_mass = 2
    params.nbi.toroidal_angle = 20.0 / 180 * pi

    return params
end