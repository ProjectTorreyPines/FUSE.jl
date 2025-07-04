"""
    case_parameters(:DTT)

DTT
"""
function case_parameters(::Type{Val{:DTT}})
    ini = ParametersInits()
    act = ParametersActors()

    ini.general.casename = "DTT"
    ini.general.init_from = :scalars

    ini.build.symmetric = true
    ini.build.divertors = :lower

    ini.equilibrium.B0 = 5.85
    ini.equilibrium.ip = 5.5e6

    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.R0 = 2.19
    ini.equilibrium.Z0 = -0.04
    ini.equilibrium.ϵ = 0.32
    ini.equilibrium.κ = 1.82
    ini.equilibrium.δ = 0.5
    ini.equilibrium.ζ = 0.0
    ini.equilibrium.𝚶 = 0.1

    # explicitly set thickness of radial build layers
    # ==============
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 0.45
    layers[:OH] = 0.325

    layers[:hfs_TF] = 0.46
    layers[:hfs_low_temp_shield] = 0.02
    layers[:hfs_vacuum_vessel] = 0.115
    layers[:hfs_wall] = 0.05

    layers[:plasma] = 1.4

    layers[:lfs_wall] = 0.05
    layers[:lfs_vacuum_vessel] = 0.41
    layers[:lfs_low_temp_shield] = 0.02
    layers[:lfs_TF] = 0.585

    layers[:gap_cryostat] = 0.5
    layers[:cryostat] = 0.10
    ini.build.layers = layers
    # ==============

    ini.build.n_first_wall_conformal_layers = 1

    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:gap_cryostat].coils_inside = 6

    ini.oh.technology = :nb3sn_iter
    ini.pf_active.technology = :nb3sn_iter
    ini.tf.technology = :nb3sn_iter

    ini.tf.shape = :miller
    ini.tf.n_coils = 18

    ini.requirements.flattop_duration = 70.0

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.35
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 20E3
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Te_sep = 100.0
    ini.core_profiles.Ti_Te_ratio = 0.6
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :Ne
    ini.core_profiles.rot_core = 0.0

    ini.core_profiles.ejima = 0.4

    resize!(ini.nb_unit, 1)
    ini.nb_unit[1].power_launched = 10e6
    ini.nb_unit[1].beam_energy = 0.5e6
    ini.nb_unit[1].beam_mass = 2.0
    ini.nb_unit[1].normalized_tangency_radius = 0.9
    ini.nb_unit[1].beam_current_fraction = [1.0, 0.0, 0.0]
    ini.nb_unit[1].current_direction = :co
    ini.nb_unit[1].offaxis = false

    resize!(ini.ec_launcher, 2)
    resize!(act.ActorSimpleEC.actuator, 2)
    ini.ec_launcher[1].power_launched = 12e6 #of 32 installed
    act.ActorSimpleEC.actuator[1].rho_0 = 0.5
    ini.ec_launcher[2].power_launched = 12e6 #of 32 installed
    act.ActorSimpleEC.actuator[2].rho_0 = 0.3

    resize!(ini.ic_antenna, 1)
    ini.ic_antenna[1].power_launched = 6e6   #of 8 installed

    #### ACT ####
    act.ActorPFdesign.symmetric = true
    act.ActorPFactive.x_points_weight = 0.025

    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    return ini, act
end

function TraceCAD(::Type{Val{:DTT}})
    x_length = 5.6
    x_offset = -0.46
    y_offset = 0.18
    return TraceCAD(:DTT, x_length, x_offset, y_offset)
end
