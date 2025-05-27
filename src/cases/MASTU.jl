"""
    case_parameters(:MASTU; init_from::Symbol)

"""
function case_parameters(::Type{Val{:MASTU}}; init_from::Symbol)
    ini = FUSE.ParametersInits()
    act = FUSE.ParametersActors()
    ini.general.casename = "MASTU"
    ini.general.init_from = init_from

    if init_from == :ods
        machine_ods = joinpath("__FUSE__", "sample", "MASTU_machine.json")
        ini.ods.filename = "$(machine_ods)"
        ActorPFactive.green_model = :point
        act.ActorVerticalStability.model = false # throws an error
    end

    ini.equilibrium.boundary_from = :MXH_params
    ini.equilibrium.R0 = 0.83
    ini.equilibrium.ϵ = 0.665
    ini.equilibrium.κ = 2.36
    ini.equilibrium.B0 = -0.55
    ini.equilibrium.MXH_params = [0.836,
      0.03,
      0.66,
      2.39,
     -0.01,
     -0.04,
      0.006,
      0.036,
     -0.009,
      0.38,
      0.0964,
     -0.0897,
     -0.0044]
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 0.8e6
    ini.equilibrium.xpoints = :lower
    
    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.37
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 1e3
    ini.core_profiles.Te_ped = 300.0
    ini.core_profiles.Te_sep = 35.0
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
    
    resize!(ini.nb_unit, 2)
    ini.nb_unit[1].power_launched = 1E6
    ini.nb_unit[1].beam_energy = 80e3
    ini.nb_unit[1].beam_mass = 2.0
    ini.nb_unit[1].template_beam = :mast_onaxis
    
    ini.nb_unit[2].power_launched = 1E6
    ini.nb_unit[2].beam_energy = 80e3
    ini.nb_unit[2].beam_mass = 2.0
    ini.nb_unit[2].template_beam = :mast_offaxis
    
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