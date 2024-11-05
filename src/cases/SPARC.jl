"""
    case_parameters(:SPARC; flux_matcher::Bool=false)

CFS/MIT SPARC design
"""
function case_parameters(::Type{Val{:SPARC}}; init_from::Symbol, flux_matcher::Bool=false)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()
    ini.general.casename = "SPARC"
    ini.general.init_from = init_from

    if init_from == :ods
        machine_ods = joinpath("__FUSE__", "sample", "OS_SPARC_Device_Description.json")
        ini.ods.filename = "$(machine_ods)"
        act.ActorWholeFacility.update_build = false
        act.ActorPFactive.green_model = :point
        act.ActorVerticalStability.model = false # throws an error
    end

    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.R0 = 1.845
    ini.equilibrium.ϵ = 0.287 #a/R0
    ini.equilibrium.κ = 2.0
    ini.equilibrium.δ = 0.55
    ini.equilibrium.ζ = -0.1
    ini.equilibrium.B0 = -12.2
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 8.7e6
    ini.equilibrium.pressure_core = 2.22e6
    ini.equilibrium.xpoints = :double

    # explicitly set thickness of 
    ini.build.n_first_wall_conformal_layers = 3
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 0.38
    layers[:OH] = 0.30
    layers[:hfs_TF] = 0.36
    layers[:hfs_gap_outer_wall_TF] = 0.01
    layers[:hfs_outer_wall] = 0.03
    layers[:hfs_vacuum_vessel] = 0.06
    layers[:hfs_inner_wall] = 0.04
    layers[:plasma] = 1.35
    layers[:lfs_inner_wall] = 0.25
    layers[:lfs_vacuum_vessel] = 0.12
    layers[:lfs_outer_wall] = 0.03
    layers[:lfs_gap_outer_wall_TF] = 0.4
    layers[:lfs_TF] = 0.40
    layers[:gap_cryostat] = 0.7
    layers[:cryostat] = 0.1
    ini.build.layers = layers

    ini.build.layers[:OH].coils_inside = 8
    ini.build.layers[:lfs_inner_wall].coils_inside = collect(17:20)
    ini.build.layers[:lfs_gap_outer_wall_TF].coils_inside = collect(21:22)
    ini.build.layers[:gap_cryostat].coils_inside = collect(9:16)

    ini.oh.technology = :rebco
    ini.pf_active.technology = :rebco
    ini.tf.technology = :rebco

    ini.tf.shape = :rectangle_ellipse
    ini.tf.n_coils = 18

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    ini.requirements.flattop_duration = 10.0

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.37
    ini.core_profiles.helium_fraction = 0.1
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    resize!(ini.ic_antenna, 1)
    ini.ic_antenna[1].power_launched = 11.1 * 1e6 #25 MW maximum available, P_threshold = 21 MW

    #### ACT ####
    act.ActorPFdesign.symmetric = true

    if !flux_matcher
        act.ActorCoreTransport.model = :none
    end

    act.ActorFluxMatcher.max_iterations = 50
    act.ActorFluxMatcher.verbose = true

    act.ActorTGLF.electromagnetic = false
    act.ActorTGLF.sat_rule = :sat0
    act.ActorTGLF.model = :TJLF

    return ini, act
end

function TraceCAD(::Type{Val{:SPARC}})
    x_length = 4.66
    x_offset = -0.58
    y_offset = 0.29
    return TraceCAD(:SPARC, x_length, x_offset, y_offset)
end
