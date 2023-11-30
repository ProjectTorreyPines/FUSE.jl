"""
    case_parameters(:STEP)

STEP

Arguments:
"""
function case_parameters(::Type{Val{:STEP}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ec=1)
    act = ParametersActors()

    ini.general.casename = "STEP"
    ini.general.init_from = :ods
    ini.equilibrium.boundary_from = :scalars
    ini.ods.filename = joinpath(@__DIR__, "..", "sample", "STEP_pf_active.json")
    N = 1 
    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 1.180137712438599/N,
        :OH => 0.35554973675585266/N,
        :gap_TF_OH => 0.068920461457739/N,
        :hfs_TF => 0.3074802716013343/N,
        :hfs_gap_low_temp_shield_TF => 0.06892046145773922/N,
        :hfs_gap_vacuum_vessel_low_temp_shield => 0.06892046145773945/N,
        :hfs_vacuum_vessel_wall_outer => 0.02/N,
        :hfs_vacuum_vessel => 0.13230115364434794/N,
        :hfs_vacuum_vessel_wall_inner => 0.02/N,
        :hfs_gap_high_temp_shield_vacuum_vessel => 0.06892046145773945/N,
        :hfs_high_temp_shield => 0.2067613843732179/N,
        :hfs_blanket => 0.3316196291966129/N,
        :hfs_first_wall => 0.02/N,
        :plasma => 3.1014207655982666/N,
        :lfs_first_wall => 0.02/N,
        :lfs_blanket => 0.7654012311802618/N,
        :lfs_high_temp_shield => 0.2067613843732179/N,
        :lfs_gap_high_temp_shield_vacuum_vessel => 0.3446023072886959/N,
        :lfs_vacuum_vessel_wall_inner => 0.02/N,
        :lfs_vacuum_vessel => 0.13230115364434794/N,
        :lfs_vacuum_vessel_wall_outer => 0.02/N,
        :lfs_gap_vacuum_vessel_low_temp_shield => 0.068920461457739/N,
        :lfs_gap_low_temp_shield_TF => 0.17230115364434795/N,
        :lfs_TF => 0.3074802716013334/N,
        :gap_cryostat => 1.4221989470234107,
        :cryostat => 0.21209509270084714
    )
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = false
    ini.build.divertors = :double
    ini.build.n_first_wall_conformal_layers = 2

    ini.material.wall = "Tungsten"
    ini.material.blanket = "lithium-lead"
    ini.material.shield = "Steel, Stainless 316"

    ini.equilibrium.B0 = 3.2
    ini.equilibrium.R0 = 3.6
    ini.equilibrium.ϵ = 1/1.8
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

    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 12
    ini.pf_active.technology = :HTS

    ini.tf.n_coils = 18
    # Table 2, NF 55 (2015) 053027 - STEP TF made of high-Jc Nb3Sn, all other coils from ITER-type Nb3Sn
    ini.tf.technology = :HTS

    ini.oh.n_coils = 12
    ini.oh.technology = :HTS
    act.ActorPFcoilsOpt.symmetric = true
    act.ActorEquilibrium.symmetrize = true

    ini.ec_launcher[1].power_launched = 150.e6 #  some at rho = 0.7 with a 0.2 width
                                               # some in core 

    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 1.1

    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorTGLF.user_specified_model = "sat1_em_iter"

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
