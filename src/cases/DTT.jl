"""
    case_parameters(:DTT)

DTT
"""
function case_parameters(::Type{Val{:DTT}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_nb=1, n_ec=1, n_ic=1)
    act = ParametersActors()

    ini.general.casename = "DTT"
    ini.general.init_from = :scalars

    ini.build.symmetric = true
    ini.build.divertors = :lower

    ini.equilibrium.B0 = 5.85
    ini.equilibrium.ip = 5.5e6
    ini.equilibrium.pressure_core = 0.9e6

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

    layers[:lfs_wall] = 0.255
    layers[:lfs_vacuum_vessel] = 0.21
    layers[:lfs_low_temp_shield] = 0.02
    layers[:lfs_TF] = 0.585 # 4.078

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
    ini.core_profiles.ne_value = 0.366
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.T_ratio = 0.6
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :Ne

    ini.core_profiles.ejima = 0.4

    ini.nb_unit[1].power_launched = 10e6
    ini.nb_unit[1].beam_energy = 0.5e6
    ini.nb_unit[1].toroidal_angle = 40.0 * deg
    ini.ec_launcher[1].power_launched = 29e6 #of 32 installed
    ini.ic_antenna[1].power_launched = 6e6   #of 8 installed

    #### ACT ####
    act.ActorPFdesign.symmetric = true

    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    return ini, act
end

function TraceCAD(::Type{Val{:DTT}})
    x_length = 5.6
    x_offset = -0.46
    y_offset = 0.18
    return TraceCAD(:DTT, x_length, x_offset, y_offset)
end
