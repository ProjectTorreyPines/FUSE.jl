"""
    case_parameters(:STEP)

STEP

Arguments:
"""
function case_parameters(::Type{Val{:STEP}}; init_from::Symbol=:scalars)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ec=1)
    act = ParametersActors()
    #### INI ####

    ini.general.casename = "STEP"
    ini.general.init_from = init_from
    ini.ods.filename = joinpath(@__DIR__, "..", "sample", "STEP_pf_active.json")

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 0.2984356197352587,
        :OH => 0.13477737665463296,
        :gap_TF_OH => 0.016847172081829065,
        :hfs_TF => 0.4043321299638989,
        :hfs_gap_coils => 0.0,
        :hfs_thermal_shield => 0.03369434416365813,
        :hfs_vacuum_vessel => 0.5559566787003614,
        :hfs_blanket => 0.030541516245487195,
        :hfs_first_wall => 0.02,
        :plasma => 4.380264741275571,
        :lfs_first_wall => 0.02,
        :lfs_blanket => 0.6538868832731644,
        :lfs_vacuum_vessel => 0.6064981949458499,
        :lfs_thermal_shield => 0.13477737665463252 + 0.06738868832731626,
        :lfs_gap_coils => 3.0,
        :lfs_TF => 0.4043321299638989,
        :gap_cryostat => 1.5,
        :cryostat => 0.2
    )
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = true
    ini.build.divertors = :double
    ini.build.n_first_wall_conformal_layers = 5

    ini.material.wall = "Tungsten"
    ini.material.blanket = "lithium-lead"
    ini.material.shield = "Steel, Stainless 316"

    ini.equilibrium.B0 = 3.2
    ini.equilibrium.R0 = 3.6
    ini.equilibrium.ϵ = 1 / 1.8
    ini.equilibrium.κ = 2.93
    ini.equilibrium.δ = 0.59
    ini.equilibrium.pressure_core = 964500.0
    ini.equilibrium.ip = 20.9e6
    ini.equilibrium.xpoints = :double
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.greenwald_fraction = 1.0
    ini.core_profiles.greenwald_fraction_ped = 0.375
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0 # unkown
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne
    ini.core_profiles.helium_fraction = 0.01

    ini.oh.n_coils = 8
    ini.oh.technology = :HTS

    ini.pf_active.n_coils_inside = 10
    ini.pf_active.n_coils_outside = 0
    ini.pf_active.technology = :HTS

    ini.tf.n_coils = 12
    ini.tf.technology = :HTS
    ini.tf.shape = :rectangle

    act.ActorPFcoilsOpt.symmetric = true
    #    act.ActorEquilibrium.symmetrize = true

    ini.ec_launcher[1].power_launched = 150.e6 #  some at rho = 0.7 with a 0.2 width some in core 

    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 1.1

    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorTGLF.user_specified_model = "sat1_em_iter"

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
