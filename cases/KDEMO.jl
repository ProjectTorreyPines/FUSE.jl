"""
    case_parameters(:KDEMO)

KDEMO

Arguments:
"""
function case_parameters(::Type{Val{:KDEMO}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()

    ini.general.casename = "K-DEMO"
    ini.general.init_from = :scalars

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 1.180137712438599,
        :OH => 0.35554973675585266,
        :gap_TF_OH => 0.068920461457739,
        :hfs_TF => 0.3074802716013343,
        :hfs_gap_low_temp_shield_TF => 0.06892046145773922,
        :hfs_gap_vacuum_vessel_low_temp_shield => 0.06892046145773945,
        :hfs_vacuum_vessel_wall_outer => 0.02,
        :hfs_vacuum_vessel => 0.13230115364434794,
        :hfs_vacuum_vessel_wall_inner => 0.02,
        :hfs_gap_high_temp_shield_vacuum_vessel => 0.06892046145773945,
        :hfs_high_temp_shield => 0.2067613843732179,
        :hfs_blanket => 0.3316196291966129,
        :hfs_first_wall => 0.02,
        :plasma => 3.1014207655982666,
        :lfs_first_wall => 0.02,
        :lfs_blanket => 0.7654012311802618,
        :lfs_high_temp_shield => 0.2067613843732179,
        :lfs_gap_high_temp_shield_vacuum_vessel => 0.3446023072886959,
        :lfs_vacuum_vessel_wall_inner => 0.02,
        :lfs_vacuum_vessel => 0.13230115364434794,
        :lfs_vacuum_vessel_wall_outer => 0.02,
        :lfs_gap_vacuum_vessel_low_temp_shield => 0.068920461457739,
        :lfs_gap_low_temp_shield_TF => 0.17230115364434795,
        :lfs_TF => 0.3074802716013334,
        :gap_cryostat => 1.4221989470234107,
        :cryostat => 0.21209509270084714,
    )
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = false
    ini.build.divertors = :lower
    ini.build.n_first_wall_conformal_layers = 2

    ini.material.wall = "Tungsten"
    ini.material.blanket = "lithium-lead"
    ini.material.shield = "Steel, Stainless 316"

    ini.equilibrium.B0 = 7.5
    ini.equilibrium.R0 = 6.8
    ini.equilibrium.ϵ = 0.3088235294117647
    ini.equilibrium.κ = 1.85
    ini.equilibrium.δ = 0.485
    ini.equilibrium.pressure_core = 964500.0
    ini.equilibrium.ip = 1.2e7
    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.greenwald_fraction = 1.0
    ini.core_profiles.greenwald_fraction_ped = 0.675
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne
    ini.core_profiles.helium_fraction = 0.01

    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 6
    ini.pf_active.technology = :LTS

    ini.tf.n_coils = 18
    ini.tf.technology = :LTS

    ini.oh.n_coils = 6
    ini.oh.technology = :LTS

    ini.ec_launchers.power_launched = 5.0e7

    ini.ic_antennas.power_launched = 5.0e7

    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 1.1
    ini.requirements.power_electric_net = 400e6 # as example

    act.ActorFluxMatcher.evolve_densities = Dict(
        :He4 => :match_ne_scale, :He4_fast => :fixed,
        :DT => :quasi_neutrality, :DT_fast => :fixed,
        :electrons => :flux_match, :Ne20 => :match_ne_scale)

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
