"""
    case_parameters(:ARC; flux_matcher::Bool=false)

CFS/MIT ARC design
"""
function case_parameters(::Type{Val{:ARC}}; flux_matcher::Bool=false)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()
    ini.general.casename = "ARC"
    ini.general.init_from = :scalars

    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.R0 = 3.45
    ini.equilibrium.ϵ = 0.27
    ini.equilibrium.κ = 2.00
    ini.equilibrium.δ = 0.4
    ini.equilibrium.ζ = -0.1
    ini.equilibrium.B0 = -11.5
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 9.9e6

    # explicitly set thickness of radial build layers
    ini.build.n_first_wall_conformal_layers = 1
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 0.82
    layers[:OH] = 0.3
    layers[:hfs_TF] = 0.6
    layers[:hfs_gap_blanket_TF] = 0.186
    layers[:hfs_vacuum_vessel] = 0.166
    layers[:hfs_blanket] = 0.4
    layers[:hfs_wall] = 0.02
    layers[:plasma] = 2.05
    layers[:lfs_wall] = 0.02
    layers[:lfs_blanket] = 0.95
    layers[:lfs_vacuum_vessel] = 0.166
    layers[:lfs_gap_blanket_TF] = 0.28
    layers[:lfs_TF] = 0.55
    layers[:gap_cryostat] = 1.119
    layers[:cryostat] = 0.186
    ini.build.layers = layers
    ini.build.layers[:hfs_blanket].material = :flibe
    ini.build.layers[:lfs_blanket].material = :flibe

    ini.equilibrium.xpoints = :double

    ini.build.layers[:OH].coils_inside = 4
    ini.build.layers[:gap_cryostat].coils_inside = 6

    ini.oh.technology = :rebco
    ini.pf_active.technology = :rebco
    ini.tf.technology = :rebco

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 18

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    #ini.requirements.power_electric_net = 50E6 ?
    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 1.1

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.49 * 0.75
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 20E3
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.zeff = 1.5
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne
    ini.core_profiles.helium_fraction = 0.10
    ini.core_profiles.rot_core = 0.0

    resize!(ini.ic_antenna, 1)
    ini.ic_antenna[1].power_launched = 4 * 1e6 #rf power coupled

    ini.requirements.coil_j_margin = 0.1
    ini.requirements.coil_stress_margin = 0.1

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

function TraceCAD(::Type{Val{:ARC}})
    x_length = 7.23
    x_offset = 0.57
    y_offset = 0.05
    return TraceCAD(:ARC, x_length, x_offset, y_offset)
end
