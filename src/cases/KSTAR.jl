"""
    case_parameters(:KSTAR; init_from::Symbol)

"""
function case_parameters(::Type{Val{:KSTAR}}; init_from::Symbol)
    ini = FUSE.ParametersInits()
    act = FUSE.ParametersActors()
    ini.general.casename = "KSTAR"
    ini.general.init_from = init_from

    if init_from == :ods
        machine_ods = joinpath("__FUSE__", "sample", "KSTAR_machine.json")
        ini.ods.filename = "$(machine_ods)"
        act.ActorPFactive.green_model = :quad
        act.ActorVerticalStability.model = false # throws an error
    end

    ini.equilibrium.boundary_from = :MXH_params
    ini.equilibrium.R0 = 1.8
    ini.equilibrium.ϵ = 0.255
    ini.equilibrium.κ = 1.73
    ini.equilibrium.B0 = -1.8
    ini.equilibrium.MXH_params = [  1.777209135,
    -0.07886733400000001,
     0.25893824532924203,
     1.7395094170491383,
     0.16528634190285837,
     0.07721056062510584,
    -0.10005895575484387,
    -0.015553628081613003,
     0.01920281391545313,
     0.5415263709543849,
     0.035027031854412694,
    -0.07089515265171509,
    -0.003920221375226912]
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 0.8e6
    ini.equilibrium.xpoints = :lower
    
    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.3
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 1000.
    ini.core_profiles.Te_ped = 300.0
    ini.core_profiles.Te_sep = 40.0
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.w_ped = 0.08
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C
    ini.tf.technology = :copper
    ini.oh.technology = :copper
    ini.pf_active.technology = :copper
    resize!(ini.nb_unit, 3)
    ini.nb_unit[1].power_launched = 1E6
    ini.nb_unit[1].beam_energy = 80e3
    ini.nb_unit[1].beam_mass = 2.0
    ini.nb_unit[1].template_beam = :kstar_nb1
    
    ini.nb_unit[2].power_launched = 1E6
    ini.nb_unit[2].beam_energy = 80e3
    ini.nb_unit[2].beam_mass = 2.0
    ini.nb_unit[2].template_beam = :kstar_nb2

    ini.nb_unit[3].power_launched = 1E6
    ini.nb_unit[3].beam_energy = 80e3
    ini.nb_unit[3].beam_mass = 2.0
    ini.nb_unit[3].template_beam = :kstar_nb3
    
    act.ActorEquilibrium.model =:FRESCO
    act.ActorFluxMatcher.max_iterations = 50
    act.ActorFluxMatcher.verbose = true
    act.ActorFRESCO.relax= 0.1
    act.ActorTGLF.electromagnetic = false
    act.ActorTGLF.sat_rule = :sat0
    act.ActorTGLF.model = :TJLF
    ini.time.simulation_start = 0.0


    return ini, act
end