"""
    case_parameters(:MANTA; flux_matcher::Bool=false)

MANTA (Modular Adjustable Negative-Triangularity ARC)

https://arxiv.org/abs/2405.20243

https://burningplasma.org/resources/ref/Web_Seminars/MANTA_USBPO_Webinar_Presentation.pdf
"""
function case_parameters(::Type{Val{:MANTA}}; flux_matcher::Bool=false)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ic=1)
    act = ParametersActors()

    ini.general.casename = "MANTA"
    ini.general.init_from = :scalars

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 1.3,
        :OH => 0.33,
        :hfs_TF => 0.7,
        :hfs_vacuum_vessel => 0.166,
        :hfs_blanket => 0.9,
        :hfs_first_wall => 0.02,
        :plasma => 2.4,
        :lfs_first_wall => 0.02,
        :lfs_blanket => 0.75,
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
    ini.equilibrium.pressure_core = 1.0E6
    ini.equilibrium.ip = 10.e6
    ini.equilibrium.xpoints = :double
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.plasma_mode = :L_mode
    ini.core_profiles.ne_setting = :greenwald_fraction
    act.ActorPedestal.density_match = :ne_line
    ini.core_profiles.ne_value = 0.5
    ini.core_profiles.ne_shaping = 4.0
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 2.0
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Kr
    ini.core_profiles.helium_fraction = 0.025

    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:hfs_blanket].coils_inside = 4
    ini.build.layers[:lfs_blanket].coils_inside = 4

    ini.oh.technology = :rebco
    ini.pf_active.technology = :rebco
    ini.tf.technology = :rebco

    ini.tf.shape = :racetrack
    ini.tf.n_coils = 18

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    ini.ic_antenna[1].power_launched = 40.e6

    ini.requirements.power_electric_net = 90e6
    ini.requirements.tritium_breeding_ratio = 1.15
    ini.requirements.flattop_duration = 45.0 * 60.0

    #### ACT ####
    act.ActorPedestal.model = :WPED
    act.ActorWPED.ped_to_core_fraction = 0.3

    act.ActorFluxMatcher.max_iterations = 50
    act.ActorFluxMatcher.verbose = true
    act.ActorTGLF.electromagnetic = false
    act.ActorTGLF.sat_rule = :sat0
    act.ActorTGLF.model = :TJLF
    if !flux_matcher
        act.ActorCoreTransport.model = :none
    end

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
