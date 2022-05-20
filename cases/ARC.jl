function case_parameters(::Type{Val{:ARC}})
    ini = ParametersInit()
    act = ParametersActor()
    ini.general.casename = "ARC"
    ini.general.init_from = :scalars

    ini.equilibrium.R0 = 3.45
    ini.equilibrium.ϵ = 0.27 #a/R0
    ini.equilibrium.κ = 1.80 #kappa_a = 1.60, kappa_sep = 1.80
    ini.equilibrium.δ = 0.50 #delta_sep
    ini.equilibrium.B0 = -11.5
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 9.9e6
    ini.equilibrium.βn = 1.4
    ini.equilibrium.x_point = (3.0, 1.7) #estimate
    ini.equilibrium.symmetric = true
    act.ActorCXbuild.rebuild_wall = false
    act.ActorHFSsizing.fixed_aspect_ratio = true

    # explicitly set thickness of 
    ini.build.blanket = 1.0 #estimate
    ini.build.shield = 0.0
    ini.build.vessel = 0.5
    ini.build.n_first_wall_conformal_layers = 2
    ini.material.shield = "Tungsten"
    ini.material.blanket = "FLiBe"

    ini.pf_active.n_oh_coils = 4
    ini.pf_active.n_pf_coils_inside = 6
    ini.pf_active.n_pf_coils_outside = 4
    ini.pf_active.technology = coil_technology(:HTS)

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 18
    ini.tf.technology = coil_technology(:HTS)
    ini.oh.technology = coil_technology(:HTS)
    ini.oh.flattop_duration = 1800

    ini.core_profiles.ne_ped = 7e19 #estimate (from ITER)
    ini.core_profiles.greenwald_fraction = 0.49
    ini.core_profiles.helium_fraction = 0.10 #estimate
    ini.core_profiles.T_shaping = 1.8 #estimate (from ITER)
    ini.core_profiles.w_ped = 0.04 #estimate (from ITER)
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne #estimate (from ITER)

    ini.nbi.power_launched = 0.0
    ini.nbi.beam_energy = 0.0
    ini.ec.power_launched = 0.0
    ini.ic.power_launched = 4 * 1e6 #rf power coupled

    act.ActorPFcoilsOpt.symmetric = true #note: symmetric, but not evenly spaced

    return set_new_base!(ini), set_new_base!(act)
end