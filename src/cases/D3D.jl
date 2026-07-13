

"""
    case_parameters(::Val{:D3D}, shot::Int;
        fit_profiles::Bool=true,
        EFIT_tree::String="EFIT02",
        PROFILES_tree::String="ZIPFIT01",
        CER_analysis_type::String="CERAUTO",
        EFIT_run_id::String="",
        PROFILES_run_id::String="",
        omfit_host::String=get(ENV, "FUSE_OMFIT_HOST", "somega.gat.com"),
        omfit_root::String=get(ENV, "FUSE_OMFIT_ROOT", "/fusion/projects/theory/fuse/d3d_data_fetching/OMFIT-source"),
        omas_root::String=get(ENV, "FUSE_OMAS_ROOT", "/fusion/projects/theory/fuse/d3d_data_fetching/omas"),
        time_averaging::Float64=0.05,
        rho_averaging::Float64=0.25,
        new_impurity_match_power_rad::Symbol=:none,
        use_local_cache::Bool=false,
        use_interferometer::Bool=true,
        pull_gslite_min::Bool=false,
        offaxis_beams::Vector{String}=["15"],
        counter_beams::Vector{String}=["21"]
    )

DIII-D from experimental shot

`offaxis_beams` and `counter_beams` select which neutral-beam beamlines are treated as
off-axis / counter-Ip when filling beam geometry from templates in the `pull_gslite_min` path
(matched against the source name, e.g. "15"/"21"); everything else is co-Ip. The known
shot-dependent history of the 150° and 210° beamlines overrides this configuration
(see [`d3d_beam_template`](@ref)).
"""
function case_parameters(::Val{:D3D}, shot::Int;
    fit_profiles::Bool=true,
    EFIT_tree::String="EFIT02",
    PROFILES_tree::String="ZIPFIT01",
    CER_analysis_type::String="CERAUTO",
    EFIT_run_id::String="",
    PROFILES_run_id::String="",
    omfit_host::String=get(ENV, "FUSE_OMFIT_HOST", "somega.gat.com"),
    omfit_root::String=get(ENV, "FUSE_OMFIT_ROOT", "/fusion/projects/theory/fuse/d3d_data_fetching/OMFIT-source"),
    omas_root::String=get(ENV, "FUSE_OMAS_ROOT", "/fusion/projects/theory/fuse/d3d_data_fetching/omas"),
    time_averaging::Float64=0.05,
    rho_averaging::Float64=0.25,
    new_impurity_match_power_rad::Symbol=:none,
    use_local_cache::Bool=false,
    use_interferometer::Bool=true,
    pull_gslite_min::Bool=false,
    offaxis_beams::Vector{String}=["15"],
    counter_beams::Vector{String}=["21"]
)
    ini, act = case_parameters(Val(:D3D_machine))
    ini.general.casename = "D3D $shot"
    omfit_host_original = omfit_host

    if omfit_host == "localhost"
        phash = hash((EFIT_tree, PROFILES_tree, CER_analysis_type, ENV["USER"], omfit_root, omas_root))
        filename = "D3D_$(shot)_$(phash).h5"
        local_path = joinpath(tempdir(), ENV["USER"]*"_D3D_$(shot)")
    else
        # Resolve remote username using the ssh config
        output = read(`ssh -T -G $omfit_host`, String)
        omfit_user = nothing
        for line in split(output, '\n')
            if startswith(line, "user ")
                omfit_user = strip(split(line, ' ', limit=2)[2])
                break
            end
        end
        if isnothing(omfit_user)
            throw(ErrorException("Need to add $omfit_host to ~/.ssh"))
        end
        omfit_host = "$omfit_user@$omfit_host"
        phash = hash((EFIT_tree, PROFILES_tree, CER_analysis_type, omfit_user, omfit_root, omas_root))
        remote_path = get(ENV, "FUSE_SCRATCH", "/cscratch/"*omfit_user*"/d3d_data/$shot")
        filename = "D3D_$(shot)_$(phash).h5"
        local_path = joinpath(tempdir(), "$(omfit_user)_D3D_$(shot)")
    end
    if isdir(local_path) && !use_local_cache
        rm(local_path; recursive=true)
    end
    if !isdir(local_path)
        mkdir(local_path)
    end


    # to get user EFITs use (shot, USER01) to get (shot01, EFIT)
    if contains(EFIT_tree, "USER")
        efit_shot = parse(Int, "$(shot)$(EFIT_tree[5:end])")
        EFIT_tree = "EFIT"
    else
        efit_shot = shot
    end

    # Build OMAS command string (same for both localhost and remote execution)
    # Note: OUTPUT_FILE is a placeholder that gets replaced with actual path (local or remote)
    omas_command = "python -u $(omas_root)/omas/examples/fuse_data_export.py OUTPUT_FILE d3d $shot $EFIT_tree $PROFILES_tree --CER_ANALYSIS_TYPE=$CER_analysis_type"
    if length(EFIT_run_id) > 0
        omas_command *= " --EFIT_RUN_ID $EFIT_run_id"
    end
    if length(PROFILES_run_id) > 0
        omas_command *= " --PROFILES_RUN_ID $PROFILES_run_id"
    end

    if pull_gslite_min
        omas_command *= " --PULL_GSLITE_MIN"
        omfit_block = """
        """
    else
        omfit_block = """
    python -u $(omfit_root)/omfit/omfit.py $(omfit_root)/modules/RABBIT/SCRIPTS/rabbit_input_no_gui.py "shot=$shot" "output_path='$local_path'" > /dev/null 2> /dev/null    
        """
    end
    if omfit_host == "localhost"
        setup_block = """#!/bin/bash -l
        module purge
        module load omfit
        cd $local_path
        export PYTHONPATH=$(omas_root):\$PYTHONPATH
        """

        # Substitute OUTPUT_FILE placeholder with actual path for localhost
        omas_block = replace(omas_command, "OUTPUT_FILE" => "$local_path/$filename")
        omfit_sh = joinpath(local_path, "omfit.sh")
        open(omfit_sh, "w") do io
            return write(io, setup_block*omfit_block)
        end
        omas_sh = joinpath(local_path, "omas.sh")
        open(omas_sh, "w") do io
            return write(io, setup_block*omas_block)
        end
        Base.run(`chmod +x $omfit_sh`)
        Base.run(`chmod +x $omas_sh`)
        task = @async Base.run(`$omfit_sh`)
        Base.run(`$omas_sh`)
        wait(task)
    else
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
            module load omfit

            echo "Starting parallel tasks..." >&2

            # Run both tasks in parallel
            mkdir -p $remote_path
            cd $remote_path

            export PYTHONPATH=$(omas_root):\$PYTHONPATH

            $(omfit_block)

            $(replace(omas_command, "OUTPUT_FILE" => "$remote_path/$(filename)"))

            echo "Waiting for OMFIT D3D BEAMS data fetching to complete..." >&2
            wait
            echo "Transfering data from remote" >&2
            """
        open(joinpath(local_path, "remote_slurm.sh"), "w") do io
            return write(io, remote_slurm)
        end
        # local driver script
        local_driver = """
            #!/bin/bash

            # Use rsync to create directory if it doesn't exist and copy the script
            $(ssh_command(omfit_host, "\"mkdir -p $remote_path\""))
            $(upsync_command(omfit_host, ["$(local_path)/remote_slurm.sh"], remote_path))

            # Execute script remotely
            $(ssh_command(omfit_host, "\"module load omfit; cd $remote_path && bash remote_slurm.sh\""))

            # Retrieve results using rsync
            $(downsync_command(omfit_host, pull_gslite_min ? ["$remote_path/$(filename)"] : ["$remote_path/$(filename)", "$remote_path/nbi_ods_$shot.h5", "$remote_path/beams_$shot.dat"], local_path))
        """

        open(joinpath(local_path, "local_driver.sh"), "w") do io
            return write(io, local_driver)
        end

        # run data fetching
        @info "Connecting to $omfit_host"
        @info("Remote D3D data fetching for shot $shot")
        @info("Path on $omfit_host: $remote_path")
        @info("Path on Localhost: $local_path")
        if !isfile(joinpath(local_path, filename)) || !use_local_cache
            Base.run(`bash $local_path/local_driver.sh`)
        end
    end
    # load experimental ods
    if pull_gslite_min
        # gslite-min pull: nbi geometry comes from DIII-D templates (no nbi_ods file)
        ini.ods.filename = "$(ini.ods.filename),$(joinpath(local_path,filename))"
    else
        ini.ods.filename = "$(ini.ods.filename),$(joinpath(local_path,filename)),$(joinpath(local_path,"nbi_ods_$shot.h5"))"
    end
    @info("Loading files: $(join(map(basename,split(ini.ods.filename,","))," ; "))")
    ini.general.dd = dd1 = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)

    n_eq = length(ini.general.dd.equilibrium.time_slice)
    t_eq = n_eq >= 2 ? ini.general.dd.equilibrium.time_slice[2].time : (n_eq == 1 ? ini.general.dd.equilibrium.time_slice[1].time : -Inf)
    t_cp = length(ini.general.dd.core_profiles.profiles_1d) >= 2 ? ini.general.dd.core_profiles.profiles_1d[2].time : -Inf
    ini.time.simulation_start = max(t_eq, t_cp)

    # sanitize dd
    for nbu in dd1.nbi.unit
        # gslite-min pull provides beam power but not the energy-species power fractions
        if IMAS.hasdata(nbu.beam_power_fraction, :data)
            nbu.beam_power_fraction.data = maximum(nbu.beam_power_fraction.data; dims=2)
            nbu.beam_power_fraction.time = [0.0]
        end
    end

    if pull_gslite_min
        # apply known shot-dependent 150/210 beamline geometry on top of the user configuration
        for nbu in dd1.nbi.unit
            if isempty(nbu.beamlets_group)
                add_beam_examples!(nbu, d3d_beam_template(shot, nbu.name, offaxis_beams, counter_beams))
            end
        end

        for launcher in dd1.ec_launchers.beam
            IMAS.hasdata(launcher.power_launched, :data) && (launcher.power_launched.data = max.(launcher.power_launched.data, 0.0))
        end
        for unit in dd1.nbi.unit
            IMAS.hasdata(unit.power_launched, :data) && (unit.power_launched.data = max.(unit.power_launched.data, 0.0))
        end
    end

    if !isempty(dd1.equilibrium.time_slice)
        IMAS.flux_surfaces(dd1.equilibrium, IMAS.first_wall(dd1.wall)...)
    end

    # set time basis
    if ~pull_gslite_min
        tt = dd1.equilibrium.time
        ini.time.pulse_shedule_time_basis = range(tt[1], tt[end], 100)
    end
    # profile fitting starting from diagnostic measurements
    if fit_profiles
        n_cer = count(ch -> !isempty(ch.ion) && IMAS.hasdata(ch.ion[1].t_i, :data), dd1.charge_exchange.channel)
        if n_cer < 5 && CER_analysis_type != "CERQUICK"
            @warn "Shot $shot: only $n_cer CER channels with $CER_analysis_type — retrying with CERQUICK"
            return case_parameters(Val(:D3D), shot;
                fit_profiles, EFIT_tree, PROFILES_tree,
                CER_analysis_type="CERQUICK",
                EFIT_run_id, PROFILES_run_id,
                omfit_host=omfit_host_original, omfit_root, omas_root,
                time_averaging, rho_averaging,
                new_impurity_match_power_rad,
                use_local_cache, use_interferometer,
                pull_gslite_min, offaxis_beams, counter_beams)
        elseif n_cer < 5
            @warn "Shot $shot: only $n_cer CER channels with CERQUICK — falling back to ZIPFIT profiles"
        else
            ActorFitProfiles(dd1, act; time_averaging, rho_averaging, time_basis_ids=:equilibrium, use_interferometer)
        end
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

    # add rotation information if missing
    for cp1d in dd1.core_profiles.profiles_1d
        if ismissing(cp1d, :rotation_frequency_tor_sonic)
            cp1d.rotation_frequency_tor_sonic =
                IMAS.Hmode_profiles(0.0, ini.core_profiles.rot_core / 8, ini.core_profiles.rot_core, length(cp1d.grid.rho_tor_norm), 1.4, 1.4, 0.05)
        end
        # freeze ion rotation so it's stored as data, not recomputed via dynamic expression
        # (sonic2ωtor uses density/pressure which change during simulation, breaking replay)
        for ion in cp1d.ion
            IMAS.freeze!(ion, :rotation_frequency_tor)
        end
    end

    # add impurity to match total radiation
    if new_impurity_match_power_rad !== :none
        dd1_core_sources_old = deepcopy(dd1.core_sources)
        IMAS.new_impurity_radiation!(dd1, new_impurity_match_power_rad, dd1.summary.time, dd1.summary.global_quantities.power_radiated_inside_lcfs.value)
        dd1.core_sources = dd1_core_sources_old
    end

    # by default match line averaged density, but only when there are experimental profiles
    if !isempty(ini.general.dd.core_profiles)
        ini.core_profiles.ne_setting = :ne_line
        act.ActorPedestal.density_match = :ne_line
    end

    set_ini_act_from_ods!(ini, act)

    #### ACT ####

    for actuator in act.ActorSimpleEC.actuator
        actuator.rho_0 = missing
        actuator.ηcd_scale = 0.2 # based on comparisons with TORAY for shot 156905
        actuator.width = 0.05
    end

    return ini, act
end

"""
    d3d_beam_template(shot::Integer, name::AbstractString, offaxis_beams=["15"], counter_beams=["21"])

Map a DIII-D neutral-beam source name to a geometry template for `add_beam_examples!`,
combining the user's off-axis/counter configuration (matched against the source name,
e.g. "15"/"21") with the known shot-dependent history of the 150° and 210° beamlines.
A beam that is both off-axis and counter-Ip maps to `:d3d_counter_offaxis`.
"""
function d3d_beam_template(shot::Integer, name::AbstractString, offaxis_beams=["15"], counter_beams=["21"])
    offaxis = any(b -> occursin(b, name), offaxis_beams)
    counter = any(b -> occursin(b, name), counter_beams)

    if occursin("21", name)
        # 210° beamline: co before 123500, counter until 177300,
        # then always off-axis with user-configured direction
        if shot < 123500
            offaxis, counter = false, false
        elseif shot < 177300
            offaxis, counter = false, true
        else
            offaxis = true
        end
    elseif occursin("15", name) && shot < 143701
        # 150° beamline was on-axis before 143701
        offaxis = false
    end

    if offaxis && counter
        return :d3d_counter_offaxis
    elseif offaxis
        return :d3d_offaxis
    elseif counter
        return :d3d_counter
    else
        return :d3d_co
    end
end

"""
    case_parameters(::Val{:D3D}, ods_file::AbstractString)

DIII-D from ods file
"""
function case_parameters(::Val{:D3D}, ods_file::AbstractString)
    ini, act = case_parameters(Val(:D3D_machine))

    ini.general.casename = "D3D $ods_file"
    ini.ods.filename = "$(ini.ods.filename),$(ods_file)"

    ini.general.dd = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)
    set_ini_act_from_ods!(ini, act)

    return ini, act
end

"""
    case_parameters(::Val{:D3D}, dd::IMAS.dd)

DIII-D from dd file
"""
function case_parameters(::Val{:D3D}, dd::IMAS.dd)
    ini, act = case_parameters(Val(:D3D_machine))

    ini.general.casename = "D3D from dd"

    ini.general.dd = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)
    merge!(ini.general.dd, dd)

    IMAS.last_global_time(ini.general.dd)
    ini.time.simulation_start = dd.global_time

    set_ini_act_from_ods!(ini, act)

    return ini, act
end

"""
    case_parameters(::Val{:D3D}, scenario::Symbol)

DIII-D from sample cases
"""
function case_parameters(::Val{:D3D}, scenario::Symbol)
    filenames = Dict(
        :H_mode => "$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json")),$(joinpath("__FUSE__", "sample", "D3D_standard_Hmode.json"))",
        :L_mode => "$(joinpath("__FUSE__", "sample", "D3D_standard_Lmode.json"))",
        :default => "$(joinpath("__FUSE__", "sample", "D3D_eq_ods.json"))")

    ini, act = case_parameters(Val(:D3D), filenames[scenario])
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
        act.ActorSources.nb_model = :none
        act.ActorSources.ec_model = :none
        act.ActorSources.lh_model = :none
        act.ActorSources.ic_model = :none
        act.ActorSources.pellet_model = :none
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
    case_parameters(::Val{:D3D_machine})

Base DIII-D machine parameters that are then extended by the other `case_parameters(:D3D, ...)` functions
"""
function case_parameters(::Val{:D3D_machine})
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

function TraceCAD(::Val{:D3D})
    x_length = 3.7727
    x_offset = -0.0303
    y_offset = -0.0303
    return TraceCAD(:D3D, x_length, x_offset, y_offset)
end
