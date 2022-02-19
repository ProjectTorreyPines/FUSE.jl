struct ITERparameters end
parametersDispatcher[:ITER] = ITERparameters()

function Parameters(::ITERparameters; init_from)
    par = Parameters()
    par.general.init_from = init_from

    if init_from == :ods
        par.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "ITER_eq_ods.json")
    else
        par.equilibrium.R0 = 6.2
        par.equilibrium.ϵ = 0.32
        par.equilibrium.κ = 1.85
        par.equilibrium.δ = 0.485
        par.equilibrium.B0 = 5.3
        par.equilibrium.Z0 = 0.4
        par.equilibrium.ip = 15.E6
        par.equilibrium.βn = 2.0
        par.equilibrium.x_point = true
        par.equilibrium.symmetric = false
    end

    par.build.is_nuclear_facility = true

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

    par.nbi.beam_power = 16.7E6
    par.nbi.beam_energy = 1000e3

    return par
end