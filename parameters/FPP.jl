struct FPPparameters end
parametersDispatcher[:FPP] = FPPparameters()

function Parameters(::FPPparameters)
    params = Parameters()
    params.general.init_from = :gasc

    params.gasc.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "FPP_fBS_PBpR_scan.json")
    params.gasc.case = 59
    params.gasc.no_small_gaps = true

    params.build.is_nuclear_facility = true

    params.pf_active.n_oh_coils = 6
    params.pf_active.n_pf_coils_inside = 0
    params.pf_active.n_pf_coils_outside = 6

    params.core_profiles.rot_core = 0.0

    params.core_profiles.bulk = :DT
    params.core_profiles.impurity = :Ne

    return params
end
