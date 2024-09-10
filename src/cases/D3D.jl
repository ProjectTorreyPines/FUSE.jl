"""
    case_parameters(:D3D)

DIII-D

Arguments:

  - scenario: `:default` (loads an experimental d3d case), `:H_mode` (updates core_profiles and core_sources), `:L_mode` (updates core_profiles and core_sources)
"""
function case_parameters(::Type{Val{:D3D}}; scenario=:default, use_ods_sources=false)::Tuple{ParametersAllInits,ParametersAllActors}

    @assert scenario in (:default, :H_mode, :L_mode)

    ini = use_ods_sources ? ParametersInits() : ParametersInits(; n_nb=1)
    act = ParametersActors()

    ini.general.casename = "D3D $scenario scenario"
    ini.general.init_from = :ods
    ini.equilibrium.boundary_from = :ods

    machine_description = joinpath("__FUSE__", "sample", "D3D_machine.json")
    shot_mappings = Dict(
        :H_mode => Dict(
            :nbi_power => 4.56e6,
            :filename => "$(machine_description),$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json")),$(joinpath("__FUSE__", "sample", "D3D_standard_Hmode.json"))"
        ),
        :L_mode => Dict(
            :nbi_power => 2.4e6,
            :filename => "$(machine_description),$(joinpath("__FUSE__", "sample", "D3D_standard_Lmode.json"))"
        ),
        :default => Dict(
            :nbi_power => 5.0e6,
            :filename => "$(machine_description),$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json"))")
    )

    ini.time.simulation_start = 0.0
    ini.ods.filename = shot_mappings[scenario][:filename]
    ini.general.dd = load_ods(ini)

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_plug => 1.2,
        :hfs_TF => 1.9,
        :hfs_gap_OH_coils => 1.5,
        :hfs_wall => 0.5,
        :plasma => 0.0,
        :lfs_wall => 0.5,
        :lfs_gap_OH_coils => 1.8,
        :lfs_TF => 1.1,
        :gap_world => 1.0
    )
    ini.build.layers[:hfs_wall].material = :graphite

    ini.equilibrium.xpoints = :upper

    ini.build.n_first_wall_conformal_layers = 1
    act.ActorCXbuild.rebuild_wall = false
    ini.build.divertors = :double

    ini.oh.n_coils = 10
    ini.pf_active.n_coils_inside = 8
    ini.pf_active.n_coils_outside = 0
    ini.pf_active.technology = :copper

    ini.tf.shape = :triple_arc
    ini.tf.n_coils = 24
    ini.tf.technology = :copper

    ini.oh.technology = :copper

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.75 * 0.75
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 5E3
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C

    ini.requirements.flattop_duration = 5.0

    act.ActorPFdesign.symmetric = true

    act.ActorWholeFacility.update_build = false
    act.ActorFluxMatcher.evolve_pedestal = false
    act.ActorStabilityLimits.raise_on_breach = false


    if use_ods_sources
        act.ActorHCD.nb_model = :none
        act.ActorHCD.ec_model = :none
        act.ActorHCD.lh_model = :none
        act.ActorHCD.ic_model = :none
        act.ActorHCD.pellet_model = :none
    else
        empty!(ini.general.dd.core_sources)
        ini.nb_unit[1].power_launched = shot_mappings[scenario][:nbi_power]
        ini.nb_unit[1].beam_energy = 80e3
        ini.nb_unit[1].beam_mass = 2.0
        ini.nb_unit[1].toroidal_angle = 20.0 / 180 * pi
    end

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end

function TraceCAD(::Type{Val{:D3D}})
    x_length = 3.7727
    x_offset = -0.0303
    y_offset = -0.0303
    return TraceCAD(:D3D, x_length, x_offset, y_offset)
end