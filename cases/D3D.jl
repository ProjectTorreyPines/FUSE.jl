sample_path = joinpath(@__DIR__, "..", "sample")
shot_details = Dict(
    :H_mode => Dict(:time0 => 2.7, :filename => joinpath(sample_path, "D3D_standard_Hmode.json")),
    :L_mode => Dict(:time0 => 2.0 , :filename => joinpath(sample_path, "D3D_standard_Lmode.json")),
    :default => Dict(:time0 => 0.0 , :filename => joinpath(sample_path, "D3D_eq_ods.json"))
    )
"""
    case_parameters(:D3D)

DIII-D

Arguments:
* `scenario`: `:H_mode`, `:L_mode` or :default (loads an experimental d3d case)
"""
function case_parameters(::Type{Val{:D3D}}; scenario=:H_mode)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()

    ini.general.casename = "D3D $scenario scenario"
    ini.general.init_from = :ods
    ini.equilibrium.boundary_from = :ods

    ini.ods.filename = shot_details[scenario][:filename]
    ini.time.simulation_start = shot_details[scenario][:time0]

    ini.build.blanket = 0.0
    ini.build.shield = 0.0
    ini.build.vessel = 0.0
    ini.build.n_first_wall_conformal_layers = 2
    ini.material.wall = "Carbon, Graphite (reactor grade)"
    act.ActorCXbuild.rebuild_wall = false
    ini.build.divertors = :double

    ini.oh.n_coils = 10
    ini.pf_active.n_coils_inside = 8
    ini.pf_active.n_coils_outside = 0
    ini.pf_active.technology = :copper

    ini.tf.shape = :double_ellipse
    ini.tf.n_coils = 24
    ini.tf.technology = :copper

    ini.oh.technology = :copper

    ini.core_profiles.greenwald_fraction = 0.7
    ini.core_profiles.greenwald_fraction_ped = ini.core_profiles.greenwald_fraction * 0.75
    ini.core_profiles.helium_fraction = 0.0
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 5E3
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C

    ini.nbi.power_launched = 5E6
    ini.nbi.beam_energy = 80e3
    ini.nbi.beam_mass = 2.0
    ini.nbi.toroidal_angle = 20.0 / 180 * pi

    ini.requirements.flattop_duration = 5.0

    act.ActorPFcoilsOpt.symmetric = true
    act.ActorFluxMatcher.evolve_densities = Dict(
        :D => :quasi_neutrality,
        :electrons => :flux_match,
        :C => :match_ne_scale)

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end