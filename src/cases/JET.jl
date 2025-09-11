"""
    case_parameters(:JET)

JET D-T case
"""
function case_parameters(::Val{:JET})
    ini = ParametersInits()
    act = ParametersActors()

    #### INI ####

    ini.general.casename = "JET"
    ini.general.init_from = :ods
    ini.general.dd = IMAS.json2imas("/home/shinan/shi/JET/JET_99948_STATE0_eq.json"; show_warnings=false)

    ini.equilibrium.boundary_from = :ods

    #### ACT ####

    act.ActorHCD.ec_model = :none
    act.ActorHCD.ic_model = :none
    act.ActorHCD.lh_model = :none
    act.ActorHCD.nb_model = :none
    act.ActorFluxMatcher.max_iterations = 500
    act.ActorFluxMatcher.verbose = true

    act.ActorTGLF.electromagnetic = true
    act.ActorTGLF.sat_rule = :sat0
    act.ActorTGLF.model = :TJLF

    return ini, act
end
