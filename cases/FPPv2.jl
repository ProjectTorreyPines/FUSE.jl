"""
    case_parameters(::Type{Val{:FPPv2}})::Tuple{ParametersAllInits,ParametersAllActors}
"""
function case_parameters(::Type{Val{:FPPv2}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ec=1)
    act = ParametersActors()

    #### INI ####

    ini.general.casename = "FPP_v2"
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
    ini.build.n_first_wall_conformal_layers = 2

    ini.material.wall = "Tungsten"
    ini.material.blanket = "lithium-lead"
    ini.material.shield = "Steel, Stainless 316"

    ini.equilibrium.B0 = 4.713171689711136
    ini.equilibrium.R0 = 4.824432302041749
    ini.equilibrium.ϵ = 0.2857142857142857
    ini.equilibrium.κ = t -> ramp((t - Δt / 8) / Δt, 0.1) * 0.3826 + 0.5
    ini.equilibrium.δ = t -> ramp((t - Δt / 8) / Δt, 0.1) * 0.7
    ini.equilibrium.ζ = t -> ramp((t - Δt / 8) / Δt, 0.1) * 0.05
    ini.equilibrium.𝚶 = t -> ramp((t - Δt / 8) / Δt, 0.1) * 0.2
    ini.equilibrium.pressure_core = t -> ramp((t - Δt / 8) / Δt, 0.3) * 1.0e6 + 0.2e6
    ini.equilibrium.ip = t -> ramp(t / Δt, 0.05) * 7.0e6 + ramp((t - Δt / 2) / (Δt / 2), 0.125) * 1.0E6
    ini.equilibrium.xpoints = t -> step((t - Δt / 7) / (Δt / 2)) < 0.5 ? :none : :lower
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.greenwald_fraction = 1.0
    ini.core_profiles.greenwald_fraction_ped = 0.7
    ini.core_profiles.T_ratio = 0.825
    ini.core_profiles.T_shaping = 2.5
    ini.core_profiles.n_shaping = 2.5
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Kr
    ini.core_profiles.helium_fraction = 0.04

    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 5
    ini.pf_active.technology = :LTS

    ini.tf.n_coils = 16
    ini.tf.technology = :HTS

    ini.oh.n_coils = 6
    ini.oh.technology = :HTS

    ini.ec_launcher[1].power_launched = 2.5e7
    ini.ec_launcher[1].efficiency_conversion = 0.45
    ini.ec_launcher[1].efficiency_transmission = 0.8

    ini.requirements.power_electric_net = 2.0e8
    ini.requirements.flattop_duration = 40000.0
    ini.requirements.tritium_breeding_ratio = 1.1

    Δt = 200 # change pulse duration to change rate of change of plasma dynamics
    ini.time.pulse_shedule_time_basis = range(0.0, Δt; step=Δt / 1000)
    ini.time.simulation_start = Δt / 2

    #### ACT ####

    act.ActorStabilityLimits.models = [:model_201, :model_401]

    act.ActorPFcoilsOpt.update_equilibrium = false

    act.ActorEquilibrium.model = :TEQUILA

    act.ActorTEQUILA.relax = 0.25

    # finalize
    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
