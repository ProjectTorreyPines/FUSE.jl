function Parameters(::Type{Val{:FPP}})
    par = Parameters()
    par.general.init_from = :gasc

    par.gasc.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "FPP_fBS_PBpR_scan.json")
    par.gasc.case = 59
    par.gasc.no_small_gaps = true

    par.build.is_nuclear_facility = true

    par.pf_active.n_oh_coils = 6
    par.pf_active.n_pf_coils_inside = 0
    par.pf_active.n_pf_coils_outside = 6

    par.core_profiles.rot_core = 0.0

    par.core_profiles.bulk = :DT

    return par
end
