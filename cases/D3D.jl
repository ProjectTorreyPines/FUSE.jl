"""
    case_parameters(:D3D)

DIII-D

Arguments:

  - scenario: `:H_mode`, `:L_mode` or `:default` (loads an experimental d3d case)
"""
function case_parameters(::Type{Val{:D3D}}; scenario=:default)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_nb=1)
    act = ParametersActors()

    ini.general.casename = "D3D $scenario scenario"
    ini.general.init_from = :ods
    ini.equilibrium.boundary_from = :ods

    machine_description = joinpath("__FUSE__", "sample", "D3D_machine.json")
    shot_mappings = Dict(
        :H_mode => Dict(:time0 => 2.7, :filename => "$(machine_description),$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json")),$(joinpath("__FUSE__", "sample", "D3D_standard_Hmode.json"))"),
        :L_mode => Dict(:time0 => 2.0, :filename => "$(machine_description),$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json")),$(joinpath("__FUSE__", "sample", "D3D_standard_Lmode.json"))"),
        :default => Dict(:time0 => 1.0, :filename => "$(machine_description),$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json"))")
    )

    ini.ods.filename = shot_mappings[scenario][:filename]
    ini.time.simulation_start = shot_mappings[scenario][:time0]

    #ini.build.layers = layers_meters_from_fractions(; blanket=0.0, shield=0.0, vessel=0.0, pf_inside_tf=true, pf_outside_tf=false)
    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_plug => 1.2,
        :hfs_TF => 1.5,
        :hfs_gap_OH_coils => 1.5,
        :hfs_wall => 0.5,
        :plasma => 0.0,
        :lfs_wall => 0.5,
        :lfs_gap_OH_coils => 1.9,
        :lfs_TF => 0.75,
        :gap_world => 1.0
    )
    ini.build.layers[:hfs_wall].material = :graphite

    ini.build.n_first_wall_conformal_layers = 2
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

    ini.nb_unit[1].power_launched = 5E6
    ini.nb_unit[1].beam_energy = 80e3
    ini.nb_unit[1].beam_mass = 2.0
    ini.nb_unit[1].toroidal_angle = 20.0 / 180 * pi

    ini.requirements.flattop_duration = 5.0

    act.ActorPFdesign.symmetric = true
    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorWholeFacility.update_build = false

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end