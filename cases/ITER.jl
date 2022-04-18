function case_parameters(::Type{Val{:ITER}}; init_from)
    ini = InitParameters()
    act = ActorParameters()
    ini.general.casename = "ITER_$(init_from)"
    ini.general.init_from = init_from

    ini.build.blanket = 0.0
    ini.build.shield = 0.5
    ini.build.vessel = 0.125
    ini.material.shield = "Tungsten"

    if init_from == :ods
        ini.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "ITER_eq_ods.json")
        act.CXbuildActor.rebuild_wall = false
        act.HFSsizingActor.fixed_plasma_start_radius = true
    else
        ini.equilibrium.R0 = 6.2
        ini.equilibrium.ϵ = 0.32
        ini.equilibrium.κ = 1.85
        ini.equilibrium.δ = 0.485
        ini.equilibrium.B0 = -5.3
        ini.equilibrium.Z0 = 0.4
        ini.equilibrium.ip = 15e6
        ini.equilibrium.βn = 2.0
        ini.equilibrium.x_point = true
        ini.equilibrium.symmetric = false
        act.CXbuildActor.rebuild_wall = true
        act.HFSsizingActor.fixed_plasma_start_radius = true
    end

    # explicitly set thickness of 
    ini.build.layers = layers = DataStructures.OrderedDict()
    layers[:gap_OH] = 0.80
    layers[:OH] = 1.30
    layers[:hfs_TF] = 1.10
    layers[:gap_hfs_vacuum_vessel] = 0.30
    layers[:hfs_shield] = 0.40
    layers[:hfs_wall] = 0.06
    layers[:plasma] = 4.40
    layers[:lfs_wall] = 0.17
    layers[:lfs_shield] = 0.40
    layers[:gap_lfs_vacuum_vessel] = 1.05
    layers[:lfs_TF] = 1.10
    layers[:gap_cryostat] = 2.34
    layers[:cryostat] = 0.30

    ini.pf_active.n_oh_coils = 6
    ini.pf_active.n_pf_coils_inside = 0
    ini.pf_active.n_pf_coils_outside = 6
    ini.pf_active.technology = coil_technology(:ITER, :PF)

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 18
    ini.tf.technology = coil_technology(:ITER, :TF)

    ini.oh.technology = coil_technology(:ITER, :OH)
    ini.oh.flattop_duration = 1000

    ini.core_profiles.ne_ped = 7e19
    ini.core_profiles.n_peaking = 1.5
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.w_ped = 0.04
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    ini.nbi.power_launched = 2 * 16.7e6
    ini.nbi.beam_energy = 1e6
    ini.ec.power_launched = 2 * 10e6
    ini.ic.power_launched = 24 * 1e6

    act.PFcoilsOptActor.symmetric = true

    return set_new_base!(ini), set_new_base!(act)
end
