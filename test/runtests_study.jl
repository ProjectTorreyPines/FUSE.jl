using FUSE
using Test
using Distributed
import FUSE.CSV as CSV
import FUSE.HDF5 as HDF5
using FUSE.DataFrames
using FUSE.SimulationParameters.Distributions

@testset "study" begin
    sty, act = FUSE.study_parameters(:DatabaseGenerator)
    sty.server = "localhost"
    sty.n_workers = 2
    sty.file_save_mode = :append
    sty.save_folder = mktempdir()
    sty.n_simulations = 2
    sty.release_workers_after_run = false

    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)

    ini.equilibrium.R0 = 5.0 ↔ [1.0, 10.0]
    ini.equilibrium.Z0 = 0.0 ↔ truncated(Normal(0.0, 2.0); lower=0.0)

    dist1 = truncated(Normal(2, 1.0); lower=0.8, upper=Inf)
    dist2 = truncated(Normal(7.0, 1.0); upper=Inf)
    mixed_two_dists = MixtureModel([dist1, dist2], [0.3, 0.7])
    ini.equilibrium.B0 = 2.0 ↔ mixed_two_dists

    ini.core_profiles.ne_setting = :greenwald_fraction
    ini.core_profiles.ne_value = 0.2 ↔ (0.2, 0.5)
    ini.core_profiles.zeff = 1.1 ↔ truncated(Normal(1.5, 1.0); lower=1.0, upper=Inf)

    act.ActorPedestal.density_match = :ne_line
    act.ActorSimpleEC.actuator[1].rho_0 = 0.3 ↔ [0.29, 0.81]

    # The study must be created first, inside of which the parallel_environment is set.
    study = FUSE.StudyDatabaseGenerator(sty, ini, act)

    @everywhere import FUSE
    @everywhere ProgressMeter = FUSE.ProgressMeter

    @everywhere function workflow_DatabaseGenerator(dd::FUSE.IMAS.dd, ini::FUSE.ParametersAllInits, act::FUSE.ParametersAllActors)
        FUSE.init(dd, ini, act)
        return nothing
    end

    study.workflow = workflow_DatabaseGenerator

    @info "Running study with `:separate_folder` policy... "
    study.sty.data_storage_policy = :separate_folders
    FUSE.run(study)

    @info "Running study with `:merged_hdf5` policy... "
    study.sty.data_storage_policy = :merged_hdf5
    FUSE.run(study)

    @testset "inis and acts" begin
        Ncase = 10
        ini_list = [rand(ini) for _ in 1:Ncase]
        act_list = [rand(act) for _ in 1:Ncase]

        # Set very high zeff to make some cases fail
        ini_list[3].core_profiles.zeff = 100.0
        ini_list[6].core_profiles.zeff = 100.0

        @assert length(ini_list) == length(act_list)

        if isdir(sty.save_folder)
            rm(sty.save_folder; recursive=true)
        end
        mkdir(sty.save_folder)

        study = FUSE.StudyDatabaseGenerator(sty, ini_list, act_list)
        study.workflow = workflow_DatabaseGenerator

        @info "Running study for predefined `inis` and `acts` (w/ :separate_folders policy)..."
        study.sty.data_storage_policy = :separate_folders
        FUSE.run(study)

        @info "Running study for predefined `inis` and `acts` (w/ :merged_hdf5 policy)..."
        study.sty.data_storage_policy = :merged_hdf5
        FUSE.run(study)

        @info "Checking if two run cases produce the same output..."
        save_dir = study.sty.save_folder

        all_dirs = filter(isdir, readdir(save_dir; join=true))
        db_folders = filter(x -> any(i -> contains(x, "$(i)__"), 1:length(ini_list)), all_dirs)
        sorted_dirs = sort(db_folders; by=x -> parse(Int, split(splitpath(x)[2], "__")[1]))
        @test length(db_folders) == length(ini_list)

        dds1 = typeof(IMAS.dd())[]
        inis1 = FUSE.ParametersInits[]
        acts1 = FUSE.ParametersActors[]
        for folder in sorted_dirs
            out = FUSE.load(folder)
            push!(dds1, ismissing(out.dd) ? IMAS.dd() : out.dd)
            push!(inis1, out.ini)
            push!(acts1, out.act)
        end
        df1 = CSV.read(joinpath(save_dir, "output.csv"), DataFrame)

        db_file_names = filter(x -> contains(x, "database_") && endswith(x, ".h5"), readdir(save_dir; join=true))
        @test length(db_file_names) == 1

        out = FUSE.load_database(db_file_names[1])
        dds2 = out.dds
        inis2 = out.inis
        acts2 = out.acts
        df2 = out.df

        # Comparison
        @test dds1 == dds2
        @test all(.!diff.(inis1, inis2))
        @test all(.!diff.(acts1, acts2))

        extract_files = filter(x -> contains(x, "extract_") && endswith(x, ".csv"), readdir(save_dir; join=true))
        @test length(extract_files) == 1
        df2_from_extract_csv = coalesce.(CSV.read(extract_files[1], DataFrame), NaN)

        df2_success_case = filter(x -> x.status == "success", df2)
        @test isequal(df1[!, Not(:dir, :gen)], df2_success_case[!, Not(:dir, :case, :gparent, :status)])
        @test isequal(df2, df2_from_extract_csv)
    end

    # Manual release of workers
    workers_list = Distributed.workers()
    if !isempty(workers_list)
        Distributed.rmprocs(workers_list)
        @info "released $(length(workers_list)) workers: $workers_list"
    else
        @info "no workers to release"
    end

    if isdir(study.sty.save_folder)
        rm(study.sty.save_folder; recursive=true)
    end
end