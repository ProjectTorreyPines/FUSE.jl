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
    ini.build.divertors = :double
    ini.material.wall = "Tungsten"
    ini.material.shield = "Steel, Stainless 316"

    ini.equilibrium.B0 = 5.85
    ini.equilibrium.ip = 5.5e6

    act.ActorStabilityLimits.raise_on_breach = false
    ini.equilibrium.pressure_core = 0.9e6

    ini.equilibrium.xpoints = :double
    act.ActorEquilibrium.symmetrize = true

    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.R0 = 2.19
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ϵ = 0.32
    ini.equilibrium.κ = 1.65
    ini.equilibrium.δ = 0.45
    ini.equilibrium.ζ = -0.1
    ini.equilibrium.𝚶 = 0.0

    # explicitly set thickness of radial build layers
    # ==============
    ini.build.layers = layers = OrderedCollections.OrderedDict{Symbol,Float64}()
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
    # ==============
    ini.build.n_first_wall_conformal_layers = 2

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
    ini.nbi.beam_energy = 0.5e6
    ini.ec_launchers.power_launched = 29e6 #of 32 installed
    ini.ic_antennas.power_launched = 6e6   #of 8 installed

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