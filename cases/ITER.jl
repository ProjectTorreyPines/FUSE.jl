function Parameters(::Type{Val{:ITER}}; init_from)
    par = Parameters()
    par.general.casename = "ITER_$(init_from)"
    par.general.init_from = init_from

    par.build.is_nuclear_facility = true

    if init_from == :ods
        par.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "ITER_eq_ods.json")
        par.build.symmetric = false
    else
        par.equilibrium.R0 = 6.2
        par.equilibrium.ϵ = 0.32
        par.equilibrium.κ = 1.85
        par.equilibrium.δ = 0.485
        par.equilibrium.B0 = -5.3
        par.equilibrium.Z0 = 0.4
        par.equilibrium.ip = 15e6
        par.equilibrium.βn = 2.0
        par.equilibrium.x_point = true
        par.equilibrium.symmetric = false
        par.build.symmetric = true
    end

    Dict(
    "gap_OH"                => 0.556784,
    "OH"                    => 0.36256,
    "hfs_TF"                => 0.492911,
    "gap_hfs_vacuum_vessel" => 0.1,
    "hfs_shield"            => 0.319696,
    "hfs_blanket"           => 0.870937,
    "hfs_wall"              => 0.0581266,
    "plasma"                => 3.31269,
    "lfs_wall"              => 0.0581266,
    "lfs_blanket"           => 0.870937,
    "lfs_shield"            => 0.319696,
    "gap_lfs_vacuum_vessel" => 0.1,
    "lfs_TF"                => 0.521975,
    "gap_cryostat"          => 1.67035}

    par.pf_active.n_oh_coils = 6
    par.pf_active.n_pf_coils_inside = 0
    par.pf_active.n_pf_coils_outside = 6
    par.pf_active.technology = coil_technology(:ITER, :PF)

    par.tf.shape = :princeton_D
    par.tf.n_coils = 18
    par.tf.technology = coil_technology(:ITER, :TF)

    par.oh.technology = coil_technology(:ITER, :OH)
    par.oh.flattop_duration = 1000

    par.core_profiles.ne_ped = 9e19
    par.core_profiles.n_peaking = 1.5
    par.core_profiles.T_shaping = 1.8
    par.core_profiles.w_ped = 0.08
    par.core_profiles.zeff = 2.5
    par.core_profiles.rot_core = 0.0
    par.core_profiles.bulk = :DT
    par.core_profiles.impurity = :Ne

    par.nbi.beam_power = 2 * 16.7e6
    par.nbi.beam_energy = 1e6
    par.ec.power_launched = 2 * 10e6
    par.ic.power_launched = 24 * 1e6

    par.material.shield = "Tungsten"
    par.material.blanket = "lithium-lead"

    return set_new_base!(par)
end
