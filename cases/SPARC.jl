"""
    case_parameters(:SPARC)

CFS/MIT SPARC design
"""
function case_parameters(::Type{Val{:SPARC}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()
    ini.general.casename = "SPARC"
    ini.general.init_from = :scalars

    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.R0 = 1.85
    ini.equilibrium.ϵ = 0.308 #a/R0
    ini.equilibrium.κ = 1.97 #kappa_a = 1.75 (kappa_lower-single-null = 1.65), kappa_sep = 1.97
    ini.equilibrium.δ = 0.54 #delta_sep
    ini.equilibrium.B0 = -12.2
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 8.7e6
    ini.equilibrium.pressure_core = 2.22e6
    ini.equilibrium.xpoints_number = 2
    act.ActorCXbuild.rebuild_wall = false
    act.ActorHFSsizing.fixed_aspect_ratio = true

    # explicitly set thickness of 
    ini.build.n_first_wall_conformal_layers = 3
    ini.build.layers = layers = OrderedCollections.OrderedDict()
    layers[:gap_OH] = 0.38
    layers[:OH] = 0.30
    layers[:hfs_TF] = 0.40
    layers[:gap_hfs_coils] = 0.05
    layers[:hfs_wall] = 0.05
    layers[:plasma] = 1.35
    layers[:lfs_wall] = 0.2
    layers[:gap_lfs_coils] = 0.34
    layers[:lfs_TF] = 0.60
    layers[:gap_cryostat] = 0.7
    #ini.material.shield = "Tungsten"

    ini.pf_active.n_oh_coils = 6
    ini.pf_active.n_pf_coils_inside = 6
    ini.pf_active.n_pf_coils_outside = 8
    ini.pf_active.technology = coil_technology(:HTS)

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 18 #estimate (from ARC)
    ini.tf.technology = coil_technology(:HTS)
    ini.oh.technology = coil_technology(:HTS)

    ini.requirements.flattop_duration = 10.0

    ini.core_profiles.greenwald_fraction = 0.37
    ini.core_profiles.greenwald_fraction_ped = ini.core_profiles.greenwald_fraction * 0.75
    ini.core_profiles.helium_fraction = 0.1 #estimate
    ini.core_profiles.T_shaping = 1.8 #estimate (from ITER)
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne #estimate (from ITER)

    ini.ic_antennas.power_launched = 11.1 * 1e6 #25 MW maximum available, P_threshold = 21 MW

    act.ActorPFcoilsOpt.symmetric = true
    # act.ActorEquilibrium.model = :CHEASE

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end

function TraceCAD(::Type{Val{:SPARC}})
    x_length = 4.66
    x_offset = -0.58
    y_offset = 0.29
    TraceCAD(:SPARC, x_length, x_offset, y_offset)
end