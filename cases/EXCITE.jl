"""
    case_parameters(::Type{Val{:EXCITE}})::Tuple{ParametersAllInits,ParametersAllActors}
"""
function case_parameters(::Type{Val{:EXCITE}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ec=2)
    act = ParametersActors()

    #### INI ####

    ini.general.casename = "EXCITE"
    ini.general.init_from = :scalars

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 0.253,
        :OH => 0.154,
        :gap_TF_OH => 0.01,
        :hfs_TF => 0.212,
        :hfs_gap_low_temp_shield_TF => 0.01,
        :hfs_low_temp_shield => 0.101,
        :hfs_gap_vacuum_vessel_low_temp_shield => 0.01,
        :hfs_vacuum_vessel_wall_outer => 0.05,
        :hfs_vacuum_vessel => 0.05,
        :hfs_vacuum_vessel_wall_inner => 0.05,
        :hfs_high_temp_shield => 0.001,
        :hfs_blanket => 0.004,
        :hfs_first_wall => 0.005,
        :plasma => 1.10,
        :lfs_first_wall => 0.005,
        :lfs_blanket => 0.0004,
        :lfs_high_temp_shield => 0.001,
        :lfs_vacuum_vessel_wall_inner => 0.05,
        :lfs_vacuum_vessel => 0.05,
        :lfs_vacuum_vessel_wall_outer => 0.05,
        :lfs_gap_vacuum_vessel_low_temp_shield => 0.01,
        :lfs_low_temp_shield => 0.10,
        :lfs_gap_low_temp_shield_TF => 0.4,
        :lfs_TF => 0.212,
        :gap_cryostat => 0.5,
        :cryostat => 0.05)

    ini.build.plasma_gap = 0.05
    ini.build.symmetric = true
    ini.build.divertors = :double
    ini.build.n_first_wall_conformal_layers = 2

    ini.equilibrium.B0 = 6.0
    ini.equilibrium.R0 = 1.5
    ini.equilibrium.ϵ = 0.5/1.5
    ini.equilibrium.κ = t -> ramp((t - Δt / 8) / Δt, 0.1) * 0.3826 + 0.5
    ini.equilibrium.δ = t -> ramp((t - Δt / 8) / Δt, 0.1) * 0.7
    ini.equilibrium.ζ = t -> ramp((t - Δt / 8) / Δt, 0.1) * 0.05
    ini.equilibrium.𝚶 = t -> ramp((t - Δt / 8) / Δt, 0.1) * 0.2
    ini.equilibrium.pressure_core = t -> ramp((t - Δt / 8) / Δt, 0.3) * 1.5e6 + 0.2e6
    ini.equilibrium.ip = t -> ramp(t / Δt, 0.05) * 5.0e6 + ramp((t - Δt / 2) / (Δt / 2), 0.125) * 1.0E6
    ini.equilibrium.xpoints = t -> step((t - Δt / 7) / (Δt / 2)) < 0.5 ? :none : :double
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.greenwald_fraction = 0.42
    ini.core_profiles.greenwald_fraction_ped = 0.25
    ini.core_profiles.T_ratio = 0.825
    ini.core_profiles.T_shaping = 2.5
    ini.core_profiles.n_shaping = 1.5
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :Kr
    ini.core_profiles.helium_fraction = 0.00

    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 8
    ini.pf_active.technology = :Nb3Sn

    ini.tf.n_coils = 16
    ini.tf.technology = :HTS
    ini.tf.shape = :princeton_D

    ini.oh.n_coils = 6
    ini.oh.technology = :HTS

    ini.ec_launcher[1].power_launched = 25.0e6
    ini.ec_launcher[1].efficiency_conversion = 0.45
    ini.ec_launcher[1].efficiency_transmission = 0.8

    ini.ec_launcher[2].power_launched = 25.0e6
    ini.ec_launcher[2].efficiency_conversion = 0.45
    ini.ec_launcher[2].efficiency_transmission = 0.8

    ini.requirements.power_electric_net = 0.0e0
    ini.requirements.flattop_duration = 3600.0*12
    ini.requirements.tritium_breeding_ratio = 0.0

    Δt = 200 # change pulse duration to change rate of change of plasma dynamics
    ini.time.pulse_shedule_time_basis = range(0.0, Δt; step=Δt / 1000)
    ini.time.simulation_start = Δt / 2

    #### ACT ####

    act.ActorStabilityLimits.models = [:q95_gt_2, :κ_controllability]

    act.ActorEquilibrium.model = :TEQUILA

    act.ActorTEQUILA.relax = 0.25

    # finalize
    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
