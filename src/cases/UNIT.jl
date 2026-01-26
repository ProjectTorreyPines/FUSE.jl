"""
    case_parameters(::Val{:UNIT})

UNIT circular plasma
"""
function case_parameters(::Val{:UNIT})
    ini = FUSE.ParametersInits()
    act = FUSE.ParametersActors()

    ini.general.casename = "circular"
    ini.general.init_from = :scalars

    # Equilibrium parameters
    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.B0 = 5.0
    ini.equilibrium.R0 = 1.0
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ϵ = 0.3
    ini.equilibrium.κ = 1.0
    ini.equilibrium.δ = 0.0
    ini.equilibrium.ip = 1E6
    ini.equilibrium.xpoints = :none

    # Core_profiles parameters
    ini.core_profiles.plasma_mode = :L_mode
    ini.core_profiles.ne_setting = :ne_ped
    ini.core_profiles.ne_value = 1E19
    ini.core_profiles.ne_shaping = 1.0
    ini.core_profiles.Te_shaping = 1.0
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.Te_core = 1E3
    ini.core_profiles.zeff = 1.0 # note zeff = 1 means no impurity
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.ngrid = 101
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C

    #### ACT ####

    act.ActorTEQUILA.free_boundary = false
    act.ActorPedestal.model = :none
    act.ActorCoreTransport.model = :none

    return ini, act
end
