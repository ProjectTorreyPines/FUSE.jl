function Parameters(::Type{Val{:FPP}}; init_from::Symbol)
    gasc = GASC(joinpath(dirname(abspath(@__FILE__)), "..", "sample", "FPP_fBS_PBpR_scan.json"), 59)

    par = Parameters(gasc)
    par.general.casename = "FPP_$(init_from)"
    par.general.init_from = init_from

    par.gasc.no_small_gaps = true

    if init_from == :ods
        par.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "fpp_gasc_59_step.json")
        par.build.blanket = .9
        par.build.shield = 0.5
        par.build.vessel = 0.125
    else
        par.core_profiles.rot_core = 0.0
        par.core_profiles.bulk = :DT
    end

    par.build.symmetric = false

    par.tf.shape = :triple_arc
    par.tf.n_coils = 16

    par.pf_active.n_oh_coils = 6
    par.pf_active.n_pf_coils_inside = 0
    par.pf_active.n_pf_coils_outside = 4

    par.material.shield = "Tungsten"
    par.material.blanket = "FLiBe"

    return set_new_base!(par)
end
