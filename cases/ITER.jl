function Parameters(::Type{Val{:ITER}}; init_from)
    par = Parameters()
    par.general.casename = "ITER_$(init_from)"
    par.general.init_from = init_from

    par.build.blanket = 0.0
    par.build.shield = 0.5
    par.build.vessel = 0.125
    par.material.shield = "Tungsten"

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

    # explicitly set thickness of 
    par.build.layers = layers = DataStructures.OrderedDict()
    layers[:gap_OH] = .80
    layers[:OH] = 1.30
    layers[:hfs_TF] = 1.10
    layers[:gap_hfs_vacuum_vessel] = 0.37
    layers[:hfs_shield] = 0.40
    layers[:hfs_wall] = 0.06
    layers[:plasma] = 4.51
    layers[:lfs_wall] = 0.06
    layers[:lfs_shield] = 0.40
    layers[:gap_lfs_vacuum_vessel] = 1.05
    layers[:lfs_TF] = 1.10
    layers[:gap_cryostat] = 2.34

    par.pf_active.n_oh_coils = 6
    par.pf_active.n_pf_coils_inside = 0
    par.pf_active.n_pf_coils_outside = 6
    par.pf_active.technology = coil_technology(:ITER, :PF)

    par.tf.shape = :triple_arc
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

    return set_new_base!(par)
end
