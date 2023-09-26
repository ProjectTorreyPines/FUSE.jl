"""
    case_parameters(:DTT)

DTT
"""
function case_parameters(::Type{Val{:DTT}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()
    ini.general.casename = "DTT"
    ini.general.init_from = :scalars

    ini.build.symmetric = true
    ini.material.wall = "Tungsten"
    ini.material.shield = "Steel, Stainless 316"

    ini.equilibrium.B0 = 5.85
    ini.equilibrium.ip = 5.5e6
    ini.equilibrium.pressure_core = 0.9e6 

    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = :scalars

    R0 = 2.19
    Z0 = 0.0
    ϵ = 0.32
    κ = 1.65
    δ = 0.35
    ζ = -0.09583 #lasciato quello di iter

    ini.equilibrium.R0 = R0
    ini.equilibrium.Z0 = Z0
    ini.equilibrium.ϵ = ϵ
    ini.equilibrium.κ = κ
    ini.equilibrium.δ = δ
    ini.equilibrium.ζ = ζ

    # explicitly set thickness of radial build layers
    ini.build.layers = layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 0.45 # 0.45
    layers[:OH] = 0.325  # 0.775
    layers[:hfs_TF] = 0.46 #1.235

    # layers[:hfs_vacuum_vessel] = 0.15,
    # layers[:hfs_shield] = 0.05
#roba da FPP
    layers[:hfs_gap_low_temp_shield_TF] = 0.01
    layers[:hfs_low_temp_shield] = 0.02
    layers[:hfs_gap_vacuum_vessel_low_temp_shield] = 0.01 #1.275

    # raggio vessel int = 1.275
    layers[:hfs_vacuum_vessel_wall_outer] = 0.01
    layers[:hfs_vacuum_vessel] = 0.115
    layers[:hfs_vacuum_vessel_wall_inner] = 0.01 #1.410
#roba da ITER.jl
    layers[:hfs_wall] = 0.05 #1.460
    layers[:hfs_gap_plasma_wall]= 0.03 #1.490
    layers[:plasma] = 1.4 # R0 2.19 - 2.890
    layers[:lfs_gap_wall_plasma]= 0.078 # 2.968
    layers[:lfs_wall] = 0.255 #tiene conto pure del vuoto vv/fw - 3.223

    # layers[:lfs_shield] = 0.05
    # layers[:lfs_vacuum_vessel] = 0.20

#roba da FPP
    layers[:lfs_vacuum_vessel_wall_inner] = 0.01
    layers[:lfs_vacuum_vessel] = 0.21
    layers[:lfs_vacuum_vessel_wall_outer] = 0.01 # 3.453

    layers[:lfs_gap_vacuum_vessel_low_temp_shield] = 0.01 # 3.457
    layers[:lfs_low_temp_shield] = 0.02
    layers[:lfs_gap_low_temp_shield_TF] = 0.01

    layers[:lfs_TF] = 0.585 # 4.078
    layers[:gap_cryostat] = 0.5
    layers[:cryostat] = 0.30
    ini.build.n_first_wall_conformal_layers = 1

    ini.oh.n_coils = 6
    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 6
    ini.pf_active.technology = :ITER
    act.ActorPFcoilsOpt.symmetric = true

    ini.tf.shape = :double_ellipse
    ini.tf.n_coils = 18
    ini.tf.technology = :ITER

    ini.oh.technology = :ITER
    act.ActorFluxSwing.operate_oh_at_j_crit = false

    ini.requirements.flattop_duration = 70.0

    ini.core_profiles.greenwald_fraction = 0.45
    ini.core_profiles.greenwald_fraction_ped = 0.366
    ini.core_profiles.helium_fraction = 0.0
    ini.core_profiles.T_ratio = 0.6
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :Ne

    ini.core_profiles.ejima = 0.4

    ini.nbi.power_launched = 10e6         
    ini.nbi.beam_energy = 0.5e6                   #500 keV
    ini.ec_launchers.power_launched = 29e6        #of 32 installed
    ini.ic_antennas.power_launched = 6e6          #of 8 installed

    act.ActorFluxMatcher.evolve_densities = Dict(
        :Ne => :match_ne_scale,
        :D => :quasi_neutrality,
        # :He4 => :match_ne_scale,
        # :He4_fast => :fixed,
        :electrons => :flux_match)

    consistent_ini_act!(ini, act)
    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
