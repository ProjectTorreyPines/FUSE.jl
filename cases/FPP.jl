"""
    case_parameters(:FPP; version::Symbol, init_from::Symbol)

GA 2022 FPP design

Arguments:
* `version`: `:v1` or `:v1_demount`
* `init_from`: `:scalars` or `:ods` (ODS contains equilibrium information)
"""
function case_parameters(::Type{Val{:FPP}}; version::Symbol, init_from::Symbol)::Tuple{ParametersAllInits, ParametersAllActors}
    if version == :v1
        filename = "FPPv1.0_aspectRatio3.5_PBpR35.json"
        case = 0
    elseif version == :v1_demount
        filename = "FPPv1.0_aspectRatio3.5_PBpR35_demount.json"
        case = 0
    end

    # load data from FPP GASC run
    gasc = GASC(joinpath(@__DIR__, "..", "sample", filename), case)
    ini, act = case_parameters(gasc)
    ini.general.casename = "FPP_$(version)_$(init_from)"
    ini.general.init_from = init_from

    if version == :v1
        ini.build.n_first_wall_conformal_layers = 100
    end

    if init_from == :ods
        ini.ods.filename = joinpath(@__DIR__, "..", "sample", "FPPv1.0_demount_eq.json")
        act.ActorCXbuild.rebuild_wall = true # false to use wall from ODS
        act.ActorHFSsizing.fixed_aspect_ratio = true
    end

    ini.core_profiles.bulk = :DT
    ini.core_profiles.rot_core = 0.0
    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 16

    ini.pf_active.n_oh_coils = 6
    ini.pf_active.n_pf_coils_inside = 0
    if init_from == :ods
        ini.pf_active.n_pf_coils_outside = 8
    else
        ini.pf_active.n_pf_coils_outside = 6
    end

    ini.material.shield = "Tungsten"
    ini.material.blanket = "lithium-lead"

    act.ActorPFcoilsOpt.symmetric = true

    set_new_base!(ini)
    set_new_base!(act)

    # ===========
    # Add deviations from GASC below this line
    # ===========
    
    # Changing Zeff from 1.1 to 2.0 will improve confinement significantly due to the pedestal increase!
    ini.core_profiles.zeff = 2.0
    
    # ini.equilibrium.δ *= -1 # negative triangularity

    # ini.equilibrium.ζ = 0.1 # squareness
    # act.ActorEquilibrium.model = :CHEASE

    return ini, act
end
