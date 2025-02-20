"""
    case_parameters(::Type{Val{:D3D}}, shot::Int)

DIII-D from experimental shot

NOTE: calls `python` with `import omas` package to use DIII-D IMAS mappings defined there
"""
function case_parameters(::Type{Val{:D3D}}, shot::Int; EFIT_tree::String="EFIT02", PROFILES_tree::String="ZIPFIT01")
    ini, act = case_parameters(Val{:D3D_machine})

    ini.general.casename = "D3D $shot"
    shot_ods_dir = tempdir()
    omas_py = """

    print("Importing packages")
    import time
    import omas
    from numpy import *

    tic = time.time()
    ods = omas.ODS()

    print("Fetching ec_launcher data")
    omas.omas_machine.d3d.ec_launcher_active_hardware(ods, $shot)

    print("Fetching nbi data")
    omas.omas_machine.d3d.nbi_active_hardware(ods, $shot)

    print("Fetching core_profiles data")
    omas.omas_machine.d3d.core_profiles_profile_1d(ods, $shot, PROFILES_tree="$(PROFILES_tree)")

    print("Fetching wall data")
    omas.omas_machine.d3d.wall(ods, $shot)

    print("Fetching equilibrium data")
    with ods.open('d3d', $shot, options={'EFIT_tree': '$EFIT_tree'}):
        for k in range(len(ods["equilibrium.time"])):
            ods["equilibrium.time_slice"][k]["time"]
            ods["equilibrium.time_slice"][k]["global_quantities.ip"]
            ods["equilibrium.time_slice"][k]["profiles_1d.psi"]
            ods["equilibrium.time_slice"][k]["profiles_1d.f"]
            ods["equilibrium.time_slice"][k]["profiles_1d.pressure"]
            ods["equilibrium.time_slice"][k]["profiles_2d[0].psi"]
            ods["equilibrium.time_slice"][k]["profiles_2d[0].grid.dim1"]
            ods["equilibrium.time_slice"][k]["profiles_2d[0].grid.dim2"]
            ods["equilibrium.time_slice"][k]["profiles_2d[0].grid_type.index"] = 1
            ods["equilibrium.vacuum_toroidal_field.r0"]
            ods["equilibrium.vacuum_toroidal_field.b0"]

    toc = time.time()
    print(f"Data fetched in {toc-tic} seconds")

    print("Saving ODS to json")
    tic = time.time()
    ods.save("$shot_ods_dir/D3D_$shot.json")
    toc = time.time()
    print(f"Saved in {toc-tic} seconds")
    """
    @info(omas_py)
    open(joinpath(shot_ods_dir, "omas_data_fetch.py"), "w") do io
        return write(io, omas_py)
    end
    Base.run(`python -u $(joinpath(shot_ods_dir,"omas_data_fetch.py"))`)

    # load experimental ods
    ini.ods.filename = "$(ini.ods.filename),$(joinpath(shot_ods_dir, "D3D_$shot.json"))"
    ini.general.dd = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)

    # set time basis
    tt = ini.general.dd.equilibrium.time
    ini.time.pulse_shedule_time_basis = range(tt[1], tt[end], 100)

    # add flux_surfaces information to experimental dd
    IMAS.flux_surfaces(ini.general.dd.equilibrium, IMAS.first_wall(ini.general.dd.wall)...)

    # add missing rotation information
    ini.core_profiles.rot_core = 5E3
    for cp1d in ini.general.dd.core_profiles.profiles_1d
        if ismissing(cp1d, :rotation_frequency_tor_sonic)
            cp1d.rotation_frequency_tor_sonic =
                IMAS.Hmode_profiles(0.0, ini.core_profiles.rot_core / 8, ini.core_profiles.rot_core, length(cp1d.grid.rho_tor_norm), 1.4, 1.4, 0.05)
        end
    end

    ini.core_profiles.ne_setting = :ne_line
    act.ActorPedestal.density_match = :ne_line

    set_ini_act_from_ods!(ini, act)

    ini.time.simulation_start = missing # to force user selection

    #### ACT ####

    for actuator in act.ActorSimpleEC.actuator
        actuator.rho_0 = missing
    end

    return ini, act
end

"""
    case_parameters(::Type{Val{:D3D}}, scenario::AbstractString)

DIII-D from ods file
"""
function case_parameters(::Type{Val{:D3D}}, ods_file::AbstractString)
    ini, act = case_parameters(Val{:D3D_machine})

    ini.general.casename = "D3D $ods_file"
    ini.ods.filename = "$(ini.ods.filename),$(ods_file)"

    ini.general.dd = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)
    set_ini_act_from_ods!(ini, act)

    return ini, act
end

"""
    case_parameters(::Type{Val{:D3D}}, dd::IMAS.dd)

DIII-D from dd file
"""
function case_parameters(::Type{Val{:D3D}}, dd::IMAS.dd)
    ini, act = case_parameters(Val{:D3D_machine})

    ini.general.casename = "D3D from dd"

    ini.general.dd = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)
    merge!(ini.general.dd, dd)
    IMAS.last_global_time(ini.general.dd)
    ini.time.simulation_start = dd.global_time

    set_ini_act_from_ods!(ini, act)

    return ini, act
end

"""
    case_parameters(::Type{Val{:D3D}}, scenario::Symbol)

DIII-D from sample cases
"""
function case_parameters(::Type{Val{:D3D}}, scenario::Symbol)
    filenames = Dict(
        :H_mode => "$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json")),$(joinpath("__FUSE__", "sample", "D3D_standard_Hmode.json"))",
        :L_mode => "$(joinpath("__FUSE__", "sample", "D3D_standard_Lmode.json"))",
        :default => "$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json"))")

    ini, act = case_parameters(Val{:D3D}, filenames[scenario])
    ini.general.casename = "D3D $scenario"

    if isempty(ini.general.dd.core_sources)
        resize!(ini.nb_unit, 1)
        ini.nb_unit[1].power_launched = 3E6
        ini.nb_unit[1].beam_energy = 80e3
        ini.nb_unit[1].beam_mass = 2.0
        ini.nb_unit[1].toroidal_angle = 18.0 * deg
    else
        act.ActorHCD.nb_model = :none
        act.ActorHCD.ec_model = :none
        act.ActorHCD.lh_model = :none
        act.ActorHCD.ic_model = :none
        act.ActorHCD.pellet_model = :none
    end

    if isempty(ini.general.dd.core_profiles)
        ini.core_profiles.ne_setting = :greenwald_fraction_ped
        ini.core_profiles.ne_value = 0.75 * 0.75
        ini.core_profiles.ne_shaping = 0.9
        ini.core_profiles.Te_shaping = 1.8
        ini.core_profiles.Ti_Te_ratio = 1.0
        ini.core_profiles.zeff = 2.0
        ini.core_profiles.bulk = :D
        ini.core_profiles.impurity = :C
        ini.core_profiles.rot_core = 5E3
    end

    return ini, act
end

"""
    case_parameters(::Type{Val{:D3D_machine}})

Base DIII-D parameters for machine
"""
function case_parameters(::Type{Val{:D3D_machine}})
    ini = ParametersInits()
    act = ParametersActors()

    ini.ods.filename = joinpath("__FUSE__", "sample", "D3D_machine.json")
    ini.general.init_from = :ods
    ini.equilibrium.boundary_from = :ods

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_plug => 1.2,
        :hfs_TF => 1.9,
        :hfs_gap_OH_coils => 1.0,
        :hfs_gap_coils => 0.5,
        :hfs_vessel => 0.2,
        :hfs_wall => 0.3,
        :plasma => 0.0,
        :lfs_wall => 0.5,
        :lfs_vessel => 0.2,
        :lfs_gap_coils => 1.6,
        :lfs_gap_OH_coils => 0.0,
        :lfs_TF => 1.1,
        :gap_world => 1.0
    )
    ini.build.layers[:hfs_wall].material = :graphite
    ini.build.n_first_wall_conformal_layers = 2

    ini.build.divertors = :double

    ini.build.layers[:hfs_gap_OH_coils].coils_inside = 6
    ini.build.layers[:hfs_gap_coils].coils_inside = [7:10; 16:19]
    ini.build.layers[:lfs_gap_coils].coils_inside = [11:15; 20:24]

    ini.oh.technology = :copper
    ini.pf_active.technology = :copper
    ini.tf.technology = :copper

    ini.tf.n_coils = 24
    ini.tf.shape = :triple_arc

    #### ACT ####

    act.ActorPFdesign.symmetric = true

    act.ActorWholeFacility.update_build = false

    act.ActorCXbuild.rebuild_wall = false

    act.ActorFluxMatcher.evolve_pedestal = false

    act.ActorPlasmaLimits.raise_on_breach = false

    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    # average pedestal height, not peak
    act.ActorEPED.ped_factor = 0.8

    return ini, act
end

function TraceCAD(::Type{Val{:D3D}})
    x_length = 3.7727
    x_offset = -0.0303
    y_offset = -0.0303
    return TraceCAD(:D3D, x_length, x_offset, y_offset)
end
