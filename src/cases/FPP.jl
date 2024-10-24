"""
    case_parameters(:FPP; flux_matcher::Bool=false)

GA's FPP design
"""
function case_parameters(::Type{Val{:FPP}}; flux_matcher::Bool=false)::Tuple{ParametersAllInits,ParametersAllActors}
    n_ec = 6
    ini = ParametersInits(; n_ec)
    act = ParametersActors()

    #### INI ####

    ini.general.casename = "FPP"
    ini.general.init_from = :scalars

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 1.180137712438599,
        :OH => 0.35554973675585266,
        :gap_TF_OH => 0.068920461457739,
        :hfs_TF => 0.3074802716013343,
        :hfs_gap_low_temp_shield_TF => 0.06892046145773922,
        :hfs_low_temp_shield => 0.4241901854016936,
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
        :lfs_low_temp_shield => 0.4241901854016943,
        :lfs_gap_low_temp_shield_TF => 0.17230115364434795,
        :lfs_TF => 0.3074802716013334,
        :gap_cryostat => 1.4221989470234107,
        :cryostat => 0.21209509270084714)
    if false
        delete!(ini.build.layers, :hfs_vacuum_vessel_wall_inner)
        delete!(ini.build.layers, :hfs_vacuum_vessel_wall_outer)
        delete!(ini.build.layers, :lfs_vacuum_vessel_wall_inner)
        delete!(ini.build.layers, :lfs_vacuum_vessel_wall_outer)
    end
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

    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 5
    ini.pf_active.technology = :nb3sn

    ini.tf.shape = :rectangle_ellipse
    ini.tf.n_coils = 16
    ini.tf.technology = :rebco

    ini.oh.n_coils = 6
    ini.oh.technology = :rebco

    total_ec_power = 90E6
    x = range(0.1, 0.8, n_ec)
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

    act.ActorStabilityLimits.models = [:q95_gt_2, :κ_controllability]

    act.ActorFluxMatcher.max_iterations = 500
    act.ActorFluxMatcher.verbose = true
    act.ActorTGLF.electromagnetic = true
    act.ActorTGLF.sat_rule = :sat0
    act.ActorTGLF.model = :TJLF
    if !flux_matcher
        act.ActorCoreTransport.model = :none
    end

    # finalize
    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
