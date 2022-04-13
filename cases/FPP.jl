function case_parameters(::Type{Val{:FPP}}; init_from::Symbol)
    gasc = GASC(joinpath(dirname(abspath(@__FILE__)), "..", "sample", "FPP_fBS_PBpR_scan.json"), 59)

    ini = InitParameters(gasc)
    act = ActorParameters()

    ini.general.casename = "FPP_$(init_from)"
    ini.general.init_from = init_from

    ini.gasc.no_small_gaps = true

    if init_from == :ods
        ini.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "fpp_gasc_59_step.json")
        ini.build.blanket = .9
        ini.build.shield = 0.5
        ini.build.vessel = 0.125
    else
        ini.core_profiles.rot_core = 0.0
        ini.core_profiles.bulk = :DT
    end

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 16

    ini.pf_active.n_oh_coils = 6
    ini.pf_active.n_pf_coils_inside = 0
    ini.pf_active.n_pf_coils_outside = 4

    ini.material.shield = "Tungsten"
    ini.material.blanket = "FLiBe"

    act.PFcoilsOptActor.symmetric = true

    return set_new_base!(ini), set_new_base!(act)
end
