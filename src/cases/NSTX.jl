"""
    case_parameters(:NSTX; init_from::Symbol)

"""
function case_parameters(::Val{:NSTX}; init_from::Symbol)
    ini = FUSE.ParametersInits()
    act = FUSE.ParametersActors()
    ini.general.casename = "NSTX"
    ini.general.init_from = init_from

    if init_from == :ods
        machine_ods = joinpath("__FUSE__", "sample", "NSTX_machine.json")
        ini.ods.filename = "$(machine_ods)"
        act.ActorCXbuild.rebuild_wall = false
        act.ActorWholeFacility.update_build = false
        act.ActorPFactive.green_model = :quad
        act.ActorVerticalStability.model = false # throws an error
    end
    
    ini.equilibrium.boundary_from = :MXH_params
    ini.equilibrium.R0 = 0.85
    ini.equilibrium.ϵ = 0.71 #a/R0
    ini.equilibrium.κ = 2.55
    ini.equilibrium.B0 = -0.48
    ini.equilibrium.MXH_params = [   0.8420864215,
     -0.011749055000000008,
      0.6724598972767072,
      2.5283024390997917,
      0.0014790406646011162,
     -0.023162692362122348,
     -0.00785385038356652,
      0.0007884001744054791,
      0.008223125812857135,
      0.5175937886461849,
     -0.06640016514800265,
     -0.050817191188883414,
      0.049435066732047686]
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 0.8e6
    ini.equilibrium.xpoints = :lower
    
    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.6
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 1e3
    ini.core_profiles.Te_ped = 500.0
    ini.core_profiles.Te_sep = 35.0
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.w_ped = 0.12
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C
    ini.tf.technology = :copper
    ini.oh.technology = :copper
    ini.pf_active.technology = :copper
    
    
    resize!(ini.nb_unit, 1)
    ini.nb_unit[1].power_launched = 1E6
    ini.nb_unit[1].beam_energy = 80e3
    ini.nb_unit[1].beam_mass = 2.0
    ini.nb_unit[1].template_beam = :nstx
    
    act.ActorEquilibrium.model =:FRESCO
    act.ActorFluxMatcher.max_iterations = 50
    act.ActorFluxMatcher.verbose = true
    act.ActorFRESCO.relax= 0.1
    act.ActorTGLF.electromagnetic = false
    act.ActorTGLF.sat_rule = :sat0
    act.ActorTGLF.model = :TJLF
    ini.time.simulation_start = 0.4
    
    return ini, act
end