"""
    case_parameters(:MASTU; init_from::Symbol)

"""
function case_parameters(::Val{:MASTU}; init_from::Symbol)
    ini = FUSE.ParametersInits()
    act = FUSE.ParametersActors()
    ini.general.casename = "MASTU"
    ini.general.init_from = init_from

    if init_from == :ods
        machine_ods = joinpath("__FUSE__", "sample", "MASTU_machine.json")
        ini.ods.filename = "$(machine_ods)"
        act.ActorCXbuild.rebuild_wall = false
        act.ActorWholeFacility.update_build = false
        act.ActorPFactive.green_model = :quad
        act.ActorVerticalStability.model = false # throws an error
    end

    ini.equilibrium.boundary_from = :MXH_params
    ini.equilibrium.R0 = 0.996630
    ini.equilibrium.ϵ = 0.543114
    ini.equilibrium.κ = 1.938625
    ini.equilibrium.B0 = -0.561750
    ini.equilibrium.MXH_params = [ 0.836657977635524,
    0.030050802044432556,
    0.63,
    2.387934019358035,
   -0.013993514162178618,
   -0.039461808063185284,
    0.006318252601694065,
    0.03583566475663035,
   -0.009081835797973238,
    0.38023377812105613,
    0.0963984614906691,
   -0.08967694186743654,
   -0.004389730106586556
  ]
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 0.744e6
    ini.equilibrium.xpoints = :upper
    
    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.533
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 1427.59
    ini.core_profiles.Te_ped = 198.78
    ini.core_profiles.Te_sep = 34.59
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.w_ped = 0.08
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.rot_core = 57441.57
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C
    ini.tf.technology = :copper
    ini.oh.technology = :copper
    ini.pf_active.technology = :copper
    
    resize!(ini.nb_unit, 2)
    ini.nb_unit[1].power_launched = 0.95E6
    ini.nb_unit[1].beam_energy = 80e3
    ini.nb_unit[1].beam_mass = 2.0
    ini.nb_unit[1].template_beam = :mast_onaxis
    
    ini.nb_unit[2].power_launched = 0.95E6
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

"""
    case_parameters(:MASTU, ods_file::AbstractString)

MAST-U from ODS file path
"""
function case_parameters(::Val{:MASTU}, ods_file::AbstractString)
    ini, act = case_parameters(Val(:MASTU); init_from=:ods)
    
    ini.general.casename = "MASTU $ods_file"
    ini.ods.filename = "$(ini.ods.filename),$(ods_file)"
    
    ini.general.dd = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)
    set_ini_act_from_ods!(ini, act)
    
    # Set simulation time from the loaded ODS
    if !isempty(ini.general.dd.core_profiles)
        IMAS.last_global_time(ini.general.dd)
        ini.time.simulation_start = ini.general.dd.global_time
    end
    
    return ini, act
end

"""
    case_parameters(:MASTU, scenario::Symbol)

MAST-U from sample cases
"""
function case_parameters(::Val{:MASTU}, scenario::Symbol)
    filenames = Dict(
        :H_mode => "$(joinpath("__FUSE__", "sample", "MASTU_Hmode.json"))",
        :default => "$(joinpath("__FUSE__", "sample", "MASTU_Hmode.json"))"
    )
    
    if !haskey(filenames, scenario)
        error("Scenario :$scenario not available for MASTU. Available scenarios: $(keys(filenames))")
    end
    
    ini, act = case_parameters(Val(:MASTU), filenames[scenario])
    ini.general.casename = "MASTU $scenario"
    
    act.ActorHCD.nb_model = :none
    act.ActorHCD.ec_model = :none
    act.ActorHCD.lh_model = :none
    act.ActorHCD.ic_model = :none
    act.ActorHCD.pellet_model = :none
    
    return ini, act
end

"""
    case_parameters(:MASTU, ods::IMAS.dd)

MAST-U from ODS/dd object with experimental data
"""
function case_parameters(::Val{:MASTU}, dd::IMAS.dd)
    # Get base MASTU machine parameters
    ini, act = case_parameters(Val(:MASTU); init_from=:ods)
    
    ini.general.casename = "MASTU from dd"
    
    # Load and merge the provided ODS data
    ini.general.dd = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)
    merge!(ini.general.dd, dd)
    
    # Set time from the ODS
    IMAS.last_global_time(ini.general.dd)
    ini.time.simulation_start = dd.global_time
    
    # Extract parameters from the ODS (equilibrium, profiles, etc.)
    set_ini_act_from_ods!(ini, act)
    
    return ini, act
end