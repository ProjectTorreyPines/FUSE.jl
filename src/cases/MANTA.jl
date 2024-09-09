"""
    case_parameters(:MANTA)

MANTA (Modular Adjustable Negative-Triangularity ARC)

https://arxiv.org/abs/2405.20243

https://burningplasma.org/resources/ref/Web_Seminars/MANTA_USBPO_Webinar_Presentation.pdf
"""
function case_parameters(::Type{Val{:MANTA}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ic=1)
    act = ParametersActors()

    ini.general.casename = "MANTA"
    ini.general.init_from = :scalars

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 1.3,
        :OH => 0.33,
        :hfs_TF => 0.7,
        :hfs_vacuum_vessel => 0.166,
        :hfs_blanket_coils => 0.9,
        :hfs_first_wall => 0.02,
        :plasma => 2.4,
        :lfs_first_wall => 0.02,
        :lfs_blanket_coils => 0.75,
        :lfs_vacuum_vessel => 0.166,
        :lfs_TF => 0.7,
        :gap_cryostat => 1.4,
        :cryostat => 0.2
    )
    ini.build.plasma_gap = 0.1
    ini.build.symmetric = false
    ini.build.divertors = :double
    ini.build.n_first_wall_conformal_layers = 1

    ini.equilibrium.B0 = 11.0
    ini.equilibrium.R0 = 4.55
    ini.equilibrium.ϵ = 0.2637362637
    ini.equilibrium.κ = 1.4
    ini.equilibrium.δ = -0.45
    ini.equilibrium.ζ = -0.25
    ini.equilibrium.pressure_core = 2.2E6
    ini.equilibrium.ip = 10.e6
    ini.equilibrium.xpoints = :double
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.ne_setting = :ne_ped
    ini.core_profiles.ne_value = 1.45e20
    ini.core_profiles.ne_sep_to_ped_ratio = 0.66
    ini.core_profiles.T_ratio = 1.01
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Kr
    ini.core_profiles.helium_fraction = 0.025

    ini.pf_active.n_coils_inside = 6
    ini.pf_active.n_coils_outside = 0
    ini.pf_active.technology = :rebco

    ini.tf.shape = :racetrack
    ini.tf.n_coils = 18
    ini.tf.technology = :rebco

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    ini.oh.n_coils = 6
    ini.oh.technology = :rebco

    ini.ic_antenna[1].power_launched = 40.e6

    ini.requirements.power_electric_net = 90e6
    ini.requirements.tritium_breeding_ratio = 1.15
    ini.requirements.flattop_duration = 45.0 * 60.0

    act.ActorPedestal.model = :WPED
    act.ActorWPED.ped_to_core_fraction = 0.3

    #### ACT ####
    act.ActorCoreTransport.model = :none # No flux matching possible

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end

function TraceCAD(::Type{Val{:MANTA}})
    x_length = 8.3
    x_offset = -0.4
    y_offset = -0.5
    return TraceCAD(:MANTA, x_length, x_offset, y_offset)
end
