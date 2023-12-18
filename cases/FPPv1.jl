"""
    case_parameters(::Type{Val{:FPP}}; version::Symbol, init_from::Symbol, STEP::Bool=true)::Tuple{ParametersAllInits,ParametersAllActors}

GA 2022 FPP design

Arguments:

  - `version`: `:v1` or `:v1_demount`
  - `init_from`: `:scalars` or `:ods`
  - `STEP`: plasma parameters to match STEP modeling
"""
function case_parameters(::Type{Val{:FPPv1}}; version::Symbol, init_from::Symbol, STEP::Bool=true)::Tuple{ParametersAllInits,ParametersAllActors}
    if version == :v1
        filename = "FPPv1.0_aspectRatio3.5_PBpR35.json"
        case = 0
    elseif version == :v1_demount
        filename = "FPPv1.0_aspectRatio3.5_PBpR35_demount.json"
        case = 0
    end

    # load data from FPP GASC run
    gasc = GASC(joinpath(@__DIR__, "..", "sample", filename), case)
    ini, act = case_parameters(gasc; add_wall_layers=0.02)
    ini.general.casename = "FPP_$(version)_$(init_from)"
    ini.general.init_from = init_from

    if init_from == :ods
        ini.ods.filename = joinpath(@__DIR__, "..", "sample", "highbatap_fpp_8MA_adhoc_EC.json")
        ini.time.simulation_start = 1.0
        act.ActorCXbuild.rebuild_wall = true # false to use wall from ODS
        ini.equilibrium.boundary_from = :scalars
        ini.equilibrium.xpoints = :double
        act.ActorEquilibrium.model = :TEQUILA
        act.ActorHCD.ec_model = :none
        act.ActorHCD.ic_model = :none
        act.ActorHCD.lh_model = :none
        act.ActorHCD.nb_model = :none
        act.ActorWholeFacility.update_plasma = false
        STEP = true
    end
    act.ActorTEQUILA.relax = 0.25

    ini.requirements.power_electric_net = 200e6 #W
    ini.requirements.tritium_breeding_ratio = 1.1
    ini.requirements.flattop_duration = 12 * 3600.0 # s

    ini.core_profiles.bulk = :DT
    ini.core_profiles.rot_core = 0.0
    ini.tf.shape = :double_ellipse
    ini.tf.n_coils = 16

    ini.oh.n_coils = 6
    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 5

    act.ActorPFdesign.symmetric = true

    set_new_base!(ini)
    set_new_base!(act)

    # ===========
    # Deviations from original GASC run are below this line
    # ===========

    # Changing Zeff from 1.1 to 2.0 will improve confinement significantly due to the pedestal increase!
    ini.core_profiles.zeff = 2.0

    # greenwald_fraction is a powerful knob
    ini.core_profiles.greenwald_fraction = 0.9
    ini.core_profiles.greenwald_fraction_ped = 0.75

    # κ<1.0 sets elongation as fraction of maximum controllable elongation estimate
    ini.equilibrium.κ = 0.95

    # negative triangularity
    # ini.equilibrium.δ *= -1

    # squareness
    ini.equilibrium.ζ = 0.15
    act.ActorEquilibrium.symmetrize = true

    # Based on STEP
    if STEP
        # zeff
        ini.core_profiles.zeff = 2.0
        # ech aiming
        act.ActorECsimple.rho_0 = 0.6
        # lower ip
        ini.equilibrium.ip = 8.0E6
        # higher density
        ini.core_profiles.greenwald_fraction = 1.0
        # limits (default FPP exceeds βn limits and greenwald density)
        act.ActorStabilityLimits.models = [:q95_gt_2,  :κ_controllability] # :gw_density, :beta_troyon_1984
        # update initial core
        ini.equilibrium.pressure_core = 1.5E6
    else
        act.ActorStabilityLimits.models = [:q95_gt_2,  :gw_density, :κ_controllability] # :beta_troyon_1984
    end

    # set density evolution for ActorFluxMatcher
    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorTGLF.user_specified_model = "sat2_em_d3d"

    # bucking
    ini.center_stack.bucked = false
    if ini.center_stack.bucked
        gasc_buck_OH_TF!(ini.build.layers)
    end

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
