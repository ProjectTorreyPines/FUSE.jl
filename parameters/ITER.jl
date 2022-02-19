struct ITERparameters end
parametersDispatcher[:ITER] = ITERparameters()

function Parameters(::ITERparameters; init_from)
    params = Parameters()
    params.general.init_from = init_from

    if init_from == :ods
        params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "ITER_eq_ods.json")
    else
        params.equilibrium.R0 = 6.2
        params.equilibrium.ϵ = 0.32
        params.equilibrium.κ = 1.85
        params.equilibrium.δ = 0.485
        params.equilibrium.B0 = 5.3
        params.equilibrium.Z0 = 0.4
        params.equilibrium.ip = 15.E6
        params.equilibrium.βn = 2.0
        params.equilibrium.x_point = true
        params.equilibrium.symmetric = false
    end

    params.build.is_nuclear_facility = true

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

    params.nbi.beam_power = 16.7E6
    params.nbi.beam_energy = 1000e3

    return params
end