

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
        counter_beams::Vector{String}=["21"],
        data_dir::String="",
        EFIT_fallback_tree::String="",
        EFIT_gap_factor::Float64=2.0,
        EFIT_hard_gap_s::Float64=0.035,
        efit_provenance::Union{Nothing,Dict{String,Any}}=nothing
    )

DIII-D from experimental shot

`data_dir` is where the shot data lives: `D3D_<shot>_<EFIT_tree>_<PROFILES_tree>_<CER_analysis_type>.h5`
(OMAS `fuse_data_export.py` output, see [`d3d_local_filename`](@ref)) plus `nbi_ods_<shot>.h5` (RABBIT beam
geometry) when available. With `use_local_cache=true` an existing file there is loaded as is, otherwise the
data is fetched from `omfit_host` into `data_dir`. Default: a per-user temporary directory. A pre-exported
file (e.g. produced off-site through FDP) is therefore used with `data_dir=<dir>, use_local_cache=true`.
Beams without geometry get it from [`d3d_beam_template`](@ref).

`offaxis_beams` and `counter_beams` select which neutral-beam beamlines are treated as
off-axis / counter-Ip when filling beam geometry from templates in the `pull_gslite_min` path
(matched against the source name, e.g. "15"/"21"); everything else is co-Ip. The known
shot-dependent history of the 150° and 210° beamlines overrides this configuration
(see [`d3d_beam_template`](@ref)).

`EFIT_fallback_tree` (e.g. "EFIT01") fills interior gaps of `EFIT_tree` wider than `EFIT_gap_factor` times its median
spacing (or, EFIT02 from EFIT01 only, `>= EFIT_hard_gap_s`) with the fallback tree's slices; `efit_provenance`, if given,
receives the source tree of each slice (see [`d3d_efit_fallback!`](@ref)).
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
    counter_beams::Vector{String}=["21"],
    data_dir::String="",
    EFIT_fallback_tree::String="",
    EFIT_gap_factor::Float64=2.0,
    EFIT_hard_gap_s::Float64=0.035,
    efit_provenance::Union{Nothing,Dict{String,Any}}=nothing
)
    ini, act = case_parameters(Val(:D3D_machine))
    EFIT_tree_name = EFIT_tree # EFIT_tree is rewritten for USERnn trees below
    ini.general.casename = "D3D $shot"
    omfit_host_original = omfit_host

    # where the shot data lives (fetched into it, or already there)
    if isempty(data_dir)
        user = omfit_host == "localhost" ? ENV["USER"] : d3d_omfit_user(omfit_host)
        data_dir = joinpath(tempdir(), "$(user)_D3D_$(shot)")
    end
    mkpath(data_dir)
    filename = d3d_local_filename(shot, EFIT_tree, PROFILES_tree, CER_analysis_type; EFIT_run_id, PROFILES_run_id)
    data_file = joinpath(data_dir, filename)

    if use_local_cache && isfile(data_file)
        @info("D3D shot $shot: using $data_file")
    else
        local_path = data_dir
        if omfit_host != "localhost"
            omfit_user = d3d_omfit_user(omfit_host)
            omfit_host = "$omfit_user@$omfit_host"
            remote_path = get(ENV, "FUSE_SCRATCH", "/cscratch/" * omfit_user * "/d3d_data/$shot")
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
            if omfit_host == "localhost"
                omfit_block = """
            python -u $(omfit_root)/omfit/omfit.py $(omfit_root)/modules/RABBIT/SCRIPTS/rabbit_input_no_gui.py "shot=$shot" "output_path='$local_path'" > /dev/null 2> /dev/null    
                """
            else
                omfit_block = """
            python -u $(omfit_root)/omfit/omfit.py $(omfit_root)/modules/RABBIT/SCRIPTS/rabbit_input_no_gui.py "shot=$shot" "output_path='$remote_path'" > /dev/null 2> /dev/null    
                """
            end
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
            try
                Base.run(`bash $local_path/local_driver.sh`)
            catch e
                error("D3D shot $shot: data fetching from $omfit_host_original failed ($(sprint(showerror, e))). " *
                      "Off-site, pass `data_dir` holding a pre-exported `$filename` and `use_local_cache=true`.")
            end
        end
        isfile(data_file) || error("D3D shot $shot: data fetching did not produce $data_file")
    end
    # load experimental ods, plus the RABBIT beam geometry when the fetch produced it (not with pull_gslite_min,
    # nor for pre-exported data): beams without geometry get DIII-D templates below
    ini.ods.filename = "$(ini.ods.filename),$(data_file)"
    nbi_file = joinpath(data_dir, "nbi_ods_$shot.h5")
    if isfile(nbi_file)
        ini.ods.filename *= ",$(nbi_file)"
    end
    @info("Loading files: $(join(map(basename,split(ini.ods.filename,","))," ; "))")
    ini.general.dd = dd1 = load_ods(ini; error_on_missing_coordinates=false, time_from_ods=true)

    # the fetch does not fill equilibrium.code.name: record the tree, then backfill its gaps from a fallback tree
    dd1.equilibrium.code.name = EFIT_tree_name
    if !isempty(EFIT_fallback_tree) && EFIT_fallback_tree != EFIT_tree_name
        fetch_fallback() = case_parameters(Val(:D3D), shot;
            fit_profiles=false, EFIT_tree=EFIT_fallback_tree, PROFILES_tree, CER_analysis_type,
            PROFILES_run_id, omfit_host=omfit_host_original, omfit_root, omas_root,
            use_local_cache=true, use_interferometer, pull_gslite_min, offaxis_beams, counter_beams, data_dir)[1].general.dd
        prov = d3d_efit_fallback!(dd1, fetch_fallback, EFIT_tree_name, EFIT_fallback_tree; gap_factor=EFIT_gap_factor, hard_gap_s=EFIT_hard_gap_s)
    else
        prov = efit_provenance === nothing ? nothing : d3d_efit_provenance(dd1, EFIT_tree_name, EFIT_fallback_tree; gap_factor=EFIT_gap_factor, hard_gap_s=NaN)
    end
    efit_provenance === nothing || merge!(efit_provenance, prov)

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

    # beams without geometry (no nbi_ods_<shot>.h5: gslite-min pull or pre-exported data):
    # apply known shot-dependent 150/210 beamline geometry on top of the user configuration
    for nbu in dd1.nbi.unit
        if isempty(nbu.beamlets_group)
            add_beam_examples!(nbu, d3d_beam_template(shot, nbu.name, offaxis_beams, counter_beams))
        end
    end

    if pull_gslite_min
        for launcher in dd1.ec_launchers.beam
            IMAS.hasdata(launcher.power_launched, :data) && (launcher.power_launched.data = max.(launcher.power_launched.data, 0.0))
        end
        for unit in dd1.nbi.unit
            IMAS.hasdata(unit.power_launched, :data) && (unit.power_launched.data = max.(unit.power_launched.data, 0.0))
        end
    end

    # The omas D3D mapping (division added in v0.92.0 as dead code, live since v0.96.0) divides the F-coil PTDATA by turns_with_sign although
    # PTDATA already is the IMAS per-turn current (gafusion/omas#449): such exports carry F-coils 55-58x too
    # small, which starves the fuse29 actuators and starts FRESCO far from the solution. Repaired here.
    repair_pf_active_omas_turns!(dd1.pf_active)

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
        # Require at least two time samples per channel, not just the presence of data: on
        # shot 207787 CERAUTO returns 7 channels with a single slice each (t = 1.94 s), which
        # clears a bare channel count but is useless to fit. time_basis_average! then
        # convolves that single sample onto the full time basis and yields all-NaN, and the
        # NaN coordinates blow up the Delaunay triangulation inside fit2d.
        n_cer = count(
            ch -> !isempty(ch.ion) && IMAS.hasdata(ch.ion[1].t_i, :data) && length(ch.ion[1].t_i.time) >= 2,
            dd1.charge_exchange.channel)
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
                pull_gslite_min, offaxis_beams, counter_beams, data_dir,
                EFIT_fallback_tree, EFIT_gap_factor, EFIT_hard_gap_s, efit_provenance)
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
        ini.core_profiles.rot_core = 30E3
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
    d3d_efit_interior_gaps(time::AbstractVector, gap_factor::Real, hard_gap_s::Real)

Intervals `(time[i], time[i+1])`, `i >= 2`, wider than `gap_factor` times the median spacing or `>= hard_gap_s`.

The first interval is excluded: the simulation starts at `time_slice[2]`, and a slice inserted before it would move the start.
"""
function d3d_efit_interior_gaps(time::AbstractVector, gap_factor::Real, hard_gap_s::Real)
    dts = sort!(diff(time))
    dt_med = dts[(length(dts)+1)÷2]
    gaps = [(time[i], time[i+1]) for i in 2:length(time)-1 if (time[i+1] - time[i] > gap_factor * dt_med || time[i+1] - time[i] >= hard_gap_s)]
    return gaps, dt_med
end

"""
    d3d_efit_splice_checks(dd1::IMAS.DD, dd_fb::IMAS.DD)

Whether slices of `dd_fb` can be copied into `dd1`: equal `vacuum_toroidal_field.r0` and equal (sign b0, sign ip).
Grid and profile-length differences are reported but not fatal.
"""
function d3d_efit_splice_checks(dd1::IMAS.DD, dd_fb::IMAS.DD)
    checks = Pair{String,String}[]
    ts1, tsf = dd1.equilibrium.time_slice, dd_fb.equilibrium.time_slice
    if isempty(ts1) || isempty(tsf)
        push!(checks, "slices_present" => "FAIL: $(length(ts1)) vs $(length(tsf))")
        return false, checks
    end
    ok = true
    r1, rf = dd1.equilibrium.vacuum_toroidal_field.r0, dd_fb.equilibrium.vacuum_toroidal_field.r0
    if isapprox(r1, rf; rtol=1e-9)
        push!(checks, "vacuum_toroidal_field.r0" => "PASS: $r1")
    else
        push!(checks, "vacuum_toroidal_field.r0" => "FAIL: $r1 vs $rf")
        ok = false
    end
    mean(x) = sum(x) / length(x)
    sb1 = Int(sign(mean(dd1.equilibrium.vacuum_toroidal_field.b0)))
    sbf = Int(sign(mean(dd_fb.equilibrium.vacuum_toroidal_field.b0)))
    si1 = Int(sign(mean([ts.global_quantities.ip for ts in ts1])))
    sif = Int(sign(mean([ts.global_quantities.ip for ts in tsf])))
    if sb1 == sbf && si1 == sif
        push!(checks, "cocos_sign_pair" => "PASS: (sign b0, sign ip) = ($sb1, $si1)")
    else
        push!(checks, "cocos_sign_pair" => "FAIL: ($sb1, $si1) vs ($sbf, $sif)")
        ok = false
    end
    g1, gf = ts1[1].profiles_2d[1].grid, tsf[1].profiles_2d[1].grid
    for (name, a, b) in (("grid.dim1", g1.dim1, gf.dim1), ("grid.dim2", g1.dim2, gf.dim2))
        same = length(a) == length(b) && isapprox(first(a), first(b)) && isapprox(last(a), last(b))
        push!(checks, "profiles_2d[0].$name" => (same ? "PASS: n=$(length(a)) [$(first(a)), $(last(a))]" :
                                                 "WARN: n=$(length(a)) [$(first(a)), $(last(a))] vs n=$(length(b)) [$(first(b)), $(last(b))]"))
    end
    n1, nf = length(ts1[1].profiles_1d.psi), length(tsf[1].profiles_1d.psi)
    push!(checks, "profiles_1d.psi length" => (n1 == nf ? "PASS: $n1" : "WARN: $n1 vs $nf"))
    return ok, checks
end

function d3d_efit_provenance(dd1::IMAS.DD, primary::AbstractString, fallback::AbstractString;
    gap_factor::Real, hard_gap_s::Real, gaps=Tuple{Float64,Float64}[], fill_counts=Int[], checks=Pair{String,String}[],
    slice_trees=fill(String(primary), length(dd1.equilibrium.time_slice)))
    return Dict{String,Any}(
        "primary_tree" => String(primary),
        "fallback_tree" => String(fallback),
        "gap_factor" => Float64(gap_factor),
        "hard_gap_s" => Float64(hard_gap_s),
        "gaps" => [Dict("t_lo" => lo, "t_hi" => hi, "n_filled" => get(fill_counts, i, 0)) for (i, (lo, hi)) in enumerate(gaps)],
        "checks" => checks,
        "slice_times" => [ts.time for ts in dd1.equilibrium.time_slice],
        "slice_trees" => collect(String, slice_trees))
end

"""
    d3d_efit_fallback!(dd1::IMAS.DD, fetch_fallback::Function, primary::AbstractString, fallback::AbstractString; gap_factor::Real, hard_gap_s::Real)

Fill interior gaps of `dd1.equilibrium` (see [`d3d_efit_interior_gaps`](@ref); `hard_gap_s` applies to EFIT02 -> EFIT01 only)
with the slices of the `fallback` tree that fall strictly inside them. With fewer than 3 primary slices, the fallback replaces the primary.

`fetch_fallback()` returns the fallback tree's `dd` and is only called when there is something to fill.
`equilibrium.vacuum_toroidal_field.b0` is rebuilt from each tree's own values and `equilibrium.code.name` becomes `"<primary>+<fallback>"`.
Returns the provenance (source tree of each slice, gaps, splice checks) as a Dict.
"""
function d3d_efit_fallback!(dd1::IMAS.DD, fetch_fallback::Function, primary::AbstractString, fallback::AbstractString; gap_factor::Real, hard_gap_s::Real)
    eq = dd1.equilibrium
    t1 = isempty(eq.time_slice) ? Float64[] : copy(eq.time)
    too_sparse = length(t1) < 3
    gaps = Tuple{Float64,Float64}[]
    hard_gap_s = (primary == "EFIT02" && fallback == "EFIT01") ? hard_gap_s : Inf
    unmerged(; kw...) = d3d_efit_provenance(dd1, primary, fallback; gap_factor, hard_gap_s, kw...)
    if too_sparse
        @warn "$primary has only $(length(t1)) equilibrium slice(s): using $fallback instead"
    else
        gaps, dt_med = d3d_efit_interior_gaps(t1, gap_factor, hard_gap_s)
        isempty(gaps) && return unmerged()
        @warn "$primary has $(length(gaps)) gap(s) beyond $(gap_factor)x its median spacing ($(round(dt_med; digits=4)) s) or $(hard_gap_s) s: " *
              join(("$(round(lo; digits=4))-$(round(hi; digits=4)) s" for (lo, hi) in gaps), ", ") * ", filling from $fallback"
    end

    dd_fb = fetch_fallback()
    # an empty primary has nothing to be consistent with
    ok, checks = isempty(t1) ? (true, Pair{String,String}[]) : d3d_efit_splice_checks(dd1, dd_fb)
    if !ok
        @error "$primary and $fallback cannot be spliced, keeping $primary as fetched:\n" * join(("    $k  $v" for (k, v) in checks), "\n")
        return unmerged(; gaps, checks)
    end

    # per-slice b0 is an expression of the parent's b0 array: snapshot each tree's values before moving slices
    tree_of, b0_of = Dict{Float64,String}(), Dict{Float64,Float64}()
    if too_sparse
        resize!(eq.time_slice, 0)
    else
        @assert length(eq.vacuum_toroidal_field.b0) == length(t1) "$primary b0 has $(length(eq.vacuum_toroidal_field.b0)) entries for $(length(t1)) times"
        for (t, b0) in zip(t1, eq.vacuum_toroidal_field.b0)
            tree_of[t] = String(primary)
            b0_of[t] = b0
        end
    end
    b0_fb = Dict{Float64,Float64}(zip(dd_fb.equilibrium.time, dd_fb.equilibrium.vacuum_toroidal_field.b0))
    fill_counts = zeros(Int, length(gaps))
    n_filled = 0
    for ts in dd_fb.equilibrium.time_slice
        gi = too_sparse ? 0 : findfirst(g -> g[1] < ts.time < g[2], gaps)
        gi === nothing && continue
        push!(eq.time_slice, deepcopy(ts))
        tree_of[ts.time] = String(fallback)
        b0_of[ts.time] = b0_fb[ts.time]
        gi > 0 && (fill_counts[gi] += 1)
        n_filled += 1
    end
    if n_filled == 0
        @warn "$fallback has no slices inside the gaps of $primary"
        return unmerged(; gaps, checks)
    end

    sort!(eq.time_slice; by=ts -> ts.time)
    eq.time = [ts.time for ts in eq.time_slice]
    eq.vacuum_toroidal_field.b0 = [b0_of[t] for t in eq.time]
    eq.code.name = "$(primary)+$(fallback)"
    return unmerged(; gaps, checks, fill_counts, slice_trees=[tree_of[ts.time] for ts in eq.time_slice])
end

"""
    d3d_local_filename(shot::Integer, EFIT_tree, PROFILES_tree, CER_analysis_type)

Name of the OMAS shot-data file that `case_parameters(:D3D, shot; data_dir)` fetches or loads, e.g.
`D3D_200000_EFIT02_ZIPFIT01_CERAUTO.h5` (run ids, when given, are appended to the tree names)
"""
function d3d_local_filename(shot::Integer, EFIT_tree::AbstractString, PROFILES_tree::AbstractString, CER_analysis_type::AbstractString;
    EFIT_run_id::AbstractString="", PROFILES_run_id::AbstractString="")
    return "D3D_$(shot)_$(EFIT_tree)$(EFIT_run_id)_$(PROFILES_tree)$(PROFILES_run_id)_$(CER_analysis_type).h5"
end

"""
    d3d_omfit_user(omfit_host::AbstractString)

Remote user name for `omfit_host` from the ssh config (`ssh -G`)
"""
function d3d_omfit_user(omfit_host::AbstractString)
    for line in split(read(`ssh -T -G $omfit_host`, String), '\n')
        if startswith(line, "user ")
            return String(strip(split(line, ' '; limit=2)[2]))
        end
    end
    throw(ErrorException("Need to add $omfit_host to ~/.ssh/config"))
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
    case_parameters(::Val{:D3D}, dd::IMAS.DD)

DIII-D from dd file
"""
function case_parameters(::Val{:D3D}, dd::IMAS.DD)
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
        ini.core_profiles.rot_core = 60E3
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

    ini.core_profiles.rot_core = 30E3

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


"""
    repair_pf_active_omas_turns!(pf_active::IMAS.pf_active; threshold::Float64=500.0)

Undo the omas D3D mapping bug (gafusion/omas#449; dead code since v0.92.0, live since v0.96.0 of 2026-04-08) that divides the F-coil PTDATA by
`turns_with_sign`: PTDATA already is the IMAS per-turn current (a free-boundary refit of the EFIT02 boundary of
shot 200000 needs F6A = -227.7 kA-turns = -4.14 kA x 55, the undivided PTDATA), so the exported F-coils are 55-58x
too small. The convention is a property of the dd, so the decision is made once for the group of coils whose first
element has |turns_with_sign| > 1: if the largest |current| among them over their traces is below `threshold` A
(a divided DIII-D F-coil set peaks near 154 A, a correct one at 6-9 kA), every coil of the group is multiplied by
its element turns (one warning listing them). The E-coil segments are single-turn elements and are never touched.
"""
function repair_pf_active_omas_turns!(pf_active::IMAS.pf_active; threshold::Float64=500.0)
    group = Tuple{eltype(pf_active.coil),Float64}[]
    imax = 0.0
    for coil in pf_active.coil
        isempty(coil.element) && continue
        # per-element turns: an F-coil is one 55-58-turn element, an E-coil segment 6-48 single-turn elements
        turns = ismissing(coil.element[1], :turns_with_sign) ? 1.0 : coil.element[1].turns_with_sign
        abs(turns) > 1.0 || continue
        (ismissing(coil.current, :data) || isempty(coil.current.data)) && continue
        push!(group, (coil, turns))
        imax = max(imax, maximum(abs, coil.current.data))
    end
    (isempty(group) || imax >= threshold) && return pf_active
    for (coil, turns) in group
        coil.current.data .*= turns
        if !ismissing(coil.current, :data_σ) && !isempty(coil.current.data_σ)
            coil.current.data_σ .*= abs(turns)
        end
    end
    @warn "pf_active: multi-turn coil currents divided by turns_with_sign (omas D3D mapping ≥ v0.96.0, gafusion/omas#449, |I|max = $(round(imax; digits=1)) A) multiplied back to the IMAS per-turn (PTDATA) convention: $(join((c.name for (c, _) in group), ", "))"
    return pf_active
end
