"""
    case_parameters(:FPP; flux_matcher::Bool=false)

GA's FPP design
"""
function case_parameters(::Type{Val{:FPP}}; flux_matcher::Bool=false)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    n_ec = 2
    ini = ParametersInits(; n_ec)
    act = ParametersActors()

    #### INI ####

    ini.general.casename = "FPP"
    ini.general.init_from = :scalars

    gaps_thickness = 0.05
    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 1.18,
        :OH => 0.35,
        :gap_TF_OH => gaps_thickness,
        :hfs_TF => 0.3,
        :hfs_gap_low_temp_shield_TF => gaps_thickness,
        :hfs_low_temp_shield => 0.42,
        :hfs_gap_vacuum_vessel_low_temp_shield => gaps_thickness,
        :hfs_vacuum_vessel_outer => 0.02,
        :hfs_gap_water => 0.13,
        :hfs_vacuum_vessel_inner => 0.02,
        :hfs_gap_high_temp_shield_vacuum_vessel => gaps_thickness,
        :hfs_high_temp_shield => 0.2,
        :hfs_blanket => 0.33,
        :hfs_first_wall => 0.02,
        :plasma => 3.1,
        :lfs_first_wall => 0.02,
        :lfs_blanket => 0.77,
        :lfs_high_temp_shield => 0.20,
        :lfs_gap_high_temp_shield_vacuum_vessel => gaps_thickness * 5,
        :lfs_vacuum_vessel_inner => 0.02,
        :lfs_gap_water => 0.13,
        :lfs_vacuum_vessel_outer => 0.02,
        :lfs_gap_vacuum_vessel_low_temp_shield => gaps_thickness * 2,
        :lfs_low_temp_shield => 0.42,
        :lfs_gap_low_temp_shield_TF => gaps_thickness * 2,
        :lfs_TF => 0.3,
        :gap_cryostat => 1.42,
        :cryostat => 0.2)
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = false
    ini.build.divertors = :lower
    ini.build.n_first_wall_conformal_layers = 1

    ini.equilibrium.B0 = 4.7
    ini.equilibrium.R0 = 4.9
    ini.equilibrium.ϵ = 0.28
    ini.equilibrium.κ = 0.8
    ini.equilibrium.δ = 0.6
    ini.equilibrium.ζ = 0.05
    ini.equilibrium.𝚶 = 0.1
    ini.equilibrium.pressure_core = 1.e6
    ini.equilibrium.ip = 8.0e6
    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.ne_setting = :greenwald_fraction
    act.ActorPedestal.density_match = :ne_line
    ini.core_profiles.ne_value = 0.9
    ini.core_profiles.ne_shaping = 2.5
    ini.core_profiles.T_ratio = 0.825
    ini.core_profiles.T_shaping = 2.5
    ini.core_profiles.Te_sep = 250.0
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Kr
    ini.core_profiles.helium_fraction = 0.01

    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:gap_cryostat].coils_inside = 5

    ini.oh.technology = :rebco
    ini.pf_active.technology = :nb3sn
    ini.tf.technology = :rebco

    ini.tf.shape = :rectangle_ellipse
    ini.tf.n_coils = 16

    total_ec_power = 90E6
    resize!(ini.ec_launcher, 6)
    x = range(0.1, 0.8, length(ini.ec_launcher))
    for (k, rho_0) in enumerate(x)
        ini.ec_launcher[k].power_launched = total_ec_power * rho_0^2 / sum(x)
        ini.ec_launcher[k].efficiency_conversion = 0.45
        ini.ec_launcher[k].efficiency_transmission = 0.8
        ini.ec_launcher[k].rho_0 = rho_0
    end

    ini.requirements.power_electric_net = 2.0e8
    ini.requirements.flattop_duration = 36000.0
    ini.requirements.tritium_breeding_ratio = 1.1

    #### ACT ####

    act.ActorLFSsizing.maintenance = :vertical

    act.ActorStabilityLimits.models = [:q95_gt_2, :κ_controllability]

    if !flux_matcher
        act.ActorCoreTransport.model = :none
    end

    act.ActorFluxMatcher.max_iterations = 500
    act.ActorFluxMatcher.verbose = true

    act.ActorTGLF.electromagnetic = true
    act.ActorTGLF.sat_rule = :sat0
    act.ActorTGLF.model = :TJLF

    return ini, act
end
