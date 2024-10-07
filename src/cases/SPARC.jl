"""
    case_parameters(:SPARC; flux_matcher::Bool=false)

CFS/MIT SPARC design
"""
function case_parameters(::Type{Val{:SPARC}}; flux_matcher::Bool=false)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ic=1)
    act = ParametersActors()
    ini.general.casename = "SPARC"
    ini.general.init_from = :scalars

    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.R0 = 1.85
    ini.equilibrium.ϵ = 0.308 #a/R0
    ini.equilibrium.κ = 1.97
    ini.equilibrium.δ = 0.45
    ini.equilibrium.B0 = -12.2
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 8.7e6
    ini.equilibrium.pressure_core = 2.22e6
    ini.equilibrium.xpoints = :double

    # explicitly set thickness of 
    ini.build.n_first_wall_conformal_layers = 1
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 0.38
    layers[:OH] = 0.30
    layers[:hfs_TF] = 0.40
    layers[:hfs_vacuum_vessel_coils] = 0.05
    layers[:hfs_wall] = 0.05
    layers[:plasma] = 1.35
    layers[:lfs_wall] = 0.2
    layers[:lfs_vacuum_vessel_coils] = 0.4
    layers[:lfs_TF] = 0.65
    layers[:gap_cryostat] = 0.7
    ini.build.layers = layers

    ini.oh.n_coils = 6
    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 8
    ini.pf_active.technology = :rebco

    ini.tf.shape = :miller
    ini.tf.n_coils = 18 #estimate (from ARC)
    ini.tf.technology = :rebco
    ini.oh.technology = :rebco

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    ini.requirements.flattop_duration = 10.0

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.37
    ini.core_profiles.helium_fraction = 0.1 #estimate
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne #estimate (from ITER)

    ini.ic_antenna[1].power_launched = 11.1 * 1e6 #25 MW maximum available, P_threshold = 21 MW

    #### ACT ####
    act.ActorPFdesign.symmetric = true

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

function TraceCAD(::Type{Val{:SPARC}})
    x_length = 4.66
    x_offset = -0.58
    y_offset = 0.29
    return TraceCAD(:SPARC, x_length, x_offset, y_offset)
end