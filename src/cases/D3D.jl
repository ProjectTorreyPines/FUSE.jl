"""
    case_parameters(::Type{Val{:D3D}}, shot::Int)

DIII-D from experimental shot

NOTE: calls `python` with `import omas` package to use DIII-D IMAS mappings defined there.
Use `ENV['OMAS_PYTHON']` to set which python executable to use.
"""
function case_parameters(
    ::Type{Val{:D3D}},
    shot::Int;
    new_impurity_match_power_rad::Symbol=:none,
    EFIT_tree::String="EFIT02",
    PROFILES_tree::String="ZIPFIT01",
    CER_analysis_type::String="CERAUTO",
    omega_user::String=get(ENV, "OMEGA_USER", ENV["USER"]),
    omega_omfit_root::String=get(ENV, "OMEGA_OMFIT_ROOT", "/fusion/projects/theory/fuse/d3d_data_fetching/OMFIT-source"),
    omega_omas_root::String=get(ENV, "OMEGA_OMAS_ROOT", "/fusion/projects/theory/fuse/d3d_data_fetching/omas"),
    use_local_cache::Bool=false
)
    ini, act = case_parameters(Val{:D3D_machine})
    ini.general.casename = "D3D $shot"

    # variables used for data fetching
    remote_omas_root = "\$OMAS_ROOT"
    if !isempty(omega_omas_root)
        remote_omas_root = omega_omas_root
    end
    remote_omfit_root = "\$OMFIT_ROOT"
    if !isempty(omega_omfit_root)
        remote_omfit_root = omega_omfit_root
    end
    remote_host = "$(omega_user)@omega.gat.com"
    remote_path = "/cscratch/$(omega_user)/d3d_data/$shot"
    filename = "D3D_$shot.h5"
    if occursin(r"omega.*.gat.com", get(ENV, "HOSTNAME", "Unknown"))
        local_path = remote_path
    else
        local_path = joinpath(tempdir(), "$(omega_user)_D3D_$(shot)")
        if isdir(local_path) && !use_local_cache
            rm(local_path; recursive=true)
        end
        if !isdir(local_path)
            mkdir(local_path)
        end
    end

    # remote omas script
    omas_py = """
        import time
        import omas
        from omas.omas_utils import printe
        from omas.machine_mappings import d3d
        from numpy import *

        ods = omas.ODS()

        tic = time.time()
        printe("- Fetching ec_launcher data")
        d3d.ec_launcher_active_hardware(ods, $shot)

        # printe("- Fetching nbi data")
        # d3d.nbi_active_hardware(ods, $shot)

        printe("- Fetching core_profiles data")
        d3d.core_profiles_profile_1d(ods, $shot, PROFILES_tree="$(PROFILES_tree)")

        printe("- Fetching wall data")
        d3d.wall(ods, $shot)

        printe("- Fetching coils data")
        d3d.pf_active_hardware(ods, $shot)
        d3d.pf_active_coil_current_data(ods, $shot)

        printe("- Fetching flux loops data")
        d3d.magnetics_floops_data(ods, $shot)

        printe("- Fetching magnetic probes data")
        d3d.magnetics_probes_data(ods, $shot)

        printe("- Fetching Thomson scattering data")
        d3d.thomson_scattering_data(ods, $shot)

        printe("- Fetching charge exchange data")
        d3d.charge_exchange_data(ods, $shot, analysis_type="$(CER_analysis_type)")

        printe("- Fetching summary data")
        d3d.summary(ods, $shot)

        printe("- Fetching equilibrium data")
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

        printe(f"Data fetched via OMAS in {time.time()-tic:.2f} [s]")

        printe("Saving ODS to $filename", end="")
        tic = time.time()
        ods.save("$filename")
        printe(f" Done in {time.time()-tic:.2f} [s]")
        """
    open(joinpath(local_path, "omas_data_fetch.py"), "w") do io
        return write(io, omas_py)
    end

    # remote bash/slurm script
    remote_slurm = """#!/bin/bash -l
        #SBATCH --job-name=fetch_d3d_omas
        #SBATCH --partition=short
        #SBATCH --cpus-per-task=1
        #SBATCH --ntasks=2
        #SBATCH --output=$remote_path/%j.out
        #SBATCH --error=$remote_path/%j.err
        #SBATCH --wait

        # Load any required modules
        module purge
        module load omfit/unstable

        echo "Starting parallel tasks..." >&2

        # Run both tasks in parallel
        cd $remote_path
        export PYTHONPATH=$(remote_omas_root):\$PYTHONPATH

        python -u $(remote_omfit_root)/omfit/omfit.py $(remote_omfit_root)/modules/RABBIT/SCRIPTS/rabbit_input_no_gui.py "shot=$shot" "output_path='$remote_path'" > /dev/null 2> /dev/null &

        python -u omas_data_fetch.py

        echo "Waiting for OMFIT D3D BEAMS data fetching to complete..." >&2
        wait
        echo "Transfering data from remote" >&2
        """
    open(joinpath(local_path, "remote_slurm.sh"), "w") do io
        return write(io, remote_slurm)
    end

    if occursin(r"omega.*.gat.com", get(ENV, "HOSTNAME", "Unknown"))
        # local driver script
        local_driver = """
            #!/bin/bash
            module load omfit; cd $remote_path && bash remote_slurm.sh
            """
    else
        # local driver script
        local_driver = """
            #!/bin/bash

            # Use rsync to create directory if it doesn't exist and copy the script
            $(ssh_command(remote_host, "\"mkdir -p $remote_path\""))
            $(upsync_command(remote_host, ["$(local_path)/remote_slurm.sh", "$(local_path)/omas_data_fetch.py"], remote_path))

            # Execute script remotely
            $(ssh_command(remote_host, "\"module load omfit; cd $remote_path && bash remote_slurm.sh\""))

            # Retrieve results using rsync
            $(downsync_command(remote_host, ["$remote_path/$(filename)", "$remote_path/nbi_ods_$shot.h5", "$remote_path/beams_$shot.dat"], local_path))
        """
    end
    open(joinpath(local_path, "local_driver.sh"), "w") do io
        return write(io, local_driver)
    end

    # run data fetching
    @info("Remote D3D data fetching for shot $shot")
    @info("Path on OMEGA: $remote_path")
    @info("Path on Localhost: $local_path")
    if !isfile(joinpath(local_path, filename)) || !use_local_cache
        Base.run(`bash $local_path/local_driver.sh`)
    end

    # load experimental ods
    ini.ods.filename = "$(ini.ods.filename),$(joinpath(local_path,filename)),$(joinpath(local_path,"nbi_ods_$shot.h5"))"
    @info("Loading files: $(join(map(basename,split(ini.ods.filename,","))," ; "))")
    ini.general.dd = dd1 = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)

    # sanitize dd
    for nbu in dd1.nbi.unit
        nbu.beam_power_fraction.data = maximum(nbu.beam_power_fraction.data, dims=2)
        nbu.beam_power_fraction.time = [0.0]
    end

    # set time basis
    tt = dd1.equilibrium.time
    ini.time.pulse_shedule_time_basis = range(tt[1], tt[end], 100)

    # add flux_surfaces information to experimental dd
    IMAS.flux_surfaces(dd1.equilibrium, IMAS.first_wall(dd1.wall)...)

    #FUSE.ActorFitProfiles(dd1,act)

    # add rotation information if missing
    for cp1d in dd1.core_profiles.profiles_1d
        if ismissing(cp1d, :rotation_frequency_tor_sonic)
            cp1d.rotation_frequency_tor_sonic =
                IMAS.Hmode_profiles(0.0, ini.core_profiles.rot_core / 8, ini.core_profiles.rot_core, length(cp1d.grid.rho_tor_norm), 1.4, 1.4, 0.05)
        end
    end

    # add impurity to match total radiation
    if new_impurity_match_power_rad !== :none
        dd1_core_sources_old = deepcopy(dd1.core_sources)
        IMAS.new_impurity_radiation!(dd1, new_impurity_match_power_rad, dd1.summary.time, dd1.summary.global_quantities.power_radiated_inside_lcfs.value)
        dd1.core_sources = dd1_core_sources_old
    end

    # by default match line averaged density
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
        resize!(ini.nb_unit, 3)
        ini.nb_unit[1].power_launched = 1E6
        ini.nb_unit[1].beam_energy = 80e3
        ini.nb_unit[1].beam_mass = 2.0
        ini.nb_unit[1].template_beam = :d3d_co

        ini.nb_unit[2].power_launched = 1E6
        ini.nb_unit[2].beam_energy = 80e3
        ini.nb_unit[2].beam_mass = 2.0
        ini.nb_unit[2].template_beam = :d3d_counter

        ini.nb_unit[3].power_launched = 1E6
        ini.nb_unit[3].beam_energy = 80e3
        ini.nb_unit[3].beam_mass = 2.0
        ini.nb_unit[3].template_beam = :d3d_offaxis

        resize!(ini.ec_launcher, 1)
        ini.ec_launcher[1].power_launched = 3E6
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

    ini.core_profiles.rot_core = 5E3

    #### ACT ####

    act.ActorEquilibrium.model = :FRESCO

    act.ActorPFdesign.symmetric = true

    act.ActorWholeFacility.update_build = false

    act.ActorCXbuild.rebuild_wall = false

    act.ActorFluxMatcher.evolve_pedestal = false

    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    Ω = 1.0 / 1E6
    act.ActorControllerIp.P = Ω * 100.0
    act.ActorControllerIp.I = Ω * 20.0
    act.ActorControllerIp.D = 0.0

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
