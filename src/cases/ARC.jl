"""
    case_parameters(:ARC)

CFS/MIT ARC design
"""
function case_parameters(::Type{Val{:ARC}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ic=1)
    act = ParametersActors()
    ini.general.casename = "ARC"
    ini.general.init_from = :scalars

    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.R0 = 3.45
    ini.equilibrium.ϵ = 0.27 #a/R0
    ini.equilibrium.κ = 1.80 #kappa_a = 1.60, kappa_sep = 1.80
    ini.equilibrium.δ = 0.50 #delta_sep
    ini.equilibrium.B0 = -11.5
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 9.9e6
    ini.equilibrium.pressure_core = 1.45e6

    # explicitly set thickness of radial build layers
    ini.build.n_first_wall_conformal_layers = 2
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 0.82
    layers[:OH] = 0.3
    layers[:hfs_TF] = 0.55
    layers[:gap_hfs_vacuum_vessel] = 0.186
    layers[:hfs_blanket] = 0.4
    layers[:hfs_wall] = 0.186
    layers[:plasma] = 2.05
    layers[:lfs_wall] = 0.186
    layers[:lfs_blanket] = 0.95
    layers[:gap_lfs_vacuum_vessel] = 0.186
    layers[:lfs_TF] = 0.55
    layers[:gap_cryostat] = 1.119
    layers[:cryostat] = 0.186
    act.ActorCXbuild.rebuild_wall = false
    ini.build.layers = layers
    ini.build.layers[:hfs_blanket].material = OrderedCollections.OrderedDict{Symbol,Float64}(:flibe => 1.0)
    ini.build.layers[:lfs_blanket].material = OrderedCollections.OrderedDict{Symbol,Float64}(:flibe => 1.0)

    ini.equilibrium.xpoints = :double

    ini.oh.n_coils = 4
    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 4
    ini.pf_active.technology = :rebco

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 18
    ini.tf.technology = :rebco
    ini.oh.technology = :rebco

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    #ini.requirements.power_electric_net = 50E6 ?
    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 1.1

    ini.core_profiles.greenwald_fraction = 0.49
    ini.core_profiles.greenwald_fraction_ped = ini.core_profiles.greenwald_fraction * 0.75
    ini.core_profiles.helium_fraction = 0.10 #estimate
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne #estimate (from ITER)

    ini.ic_antenna[1].power_launched = 4 * 1e6 #rf power coupled

    act.ActorPFdesign.symmetric = true
    act.ActorEquilibrium.symmetrize = true
    act.ActorCXbuild.rebuild_wall = true

    act.ActorHFSsizing.j_tolerance = 0.1
    act.ActorHFSsizing.stress_tolerance = 0.1

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end

function TraceCAD(::Type{Val{:ARC}})
    x_length = 7.23
    x_offset = 0.57
    y_offset = 0.05
    return TraceCAD(:ARC, x_length, x_offset, y_offset)
end