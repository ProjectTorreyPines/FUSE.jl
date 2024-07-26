"""
    case_parameters(:FIRST)

Formosa Integrated Research Spherical Tokamak
"""
function case_parameters(::Type{Val{:FIRST}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ec=0)
    act = ParametersActors()

    ini.general.casename = "FIRST"
    ini.general.init_from = :scalars

    ini.build.layers = layers_meters_from_fractions(;
        wall=1.0,
        blanket=0.0,
        shield=0.0,
        vessel=1.0,
        lfs_multiplier=3.0,
        pf_inside_tf=false,
        pf_outside_tf=true)
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = true
    ini.build.n_first_wall_conformal_layers = 1

    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.B0 = 0.8
    ini.equilibrium.R0 = 0.45
    ini.equilibrium.ϵ = 0.32 / ini.equilibrium.R0
    ini.equilibrium.κ = 2.4
    ini.equilibrium.δ = 0.5
    ini.equilibrium.pressure_core = 1E5
    ini.equilibrium.ip = 1.8E6

    ini.equilibrium.xpoints = :double

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.675
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :Ar

    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 6
    ini.pf_active.technology = :copper

    ini.tf.shape = :double_ellipse
    ini.tf.n_coils = 18
    ini.tf.technology = :copper

    ini.oh.n_coils = 6
    ini.oh.technology = :copper

    ini.requirements.flattop_duration = 2.0

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
