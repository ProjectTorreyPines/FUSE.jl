using FUSE
using Test
using Distributed

@testset "study" begin
    sty, act = FUSE.study_parameters(:DatabaseGenerator)
    sty.server = "localhost"
    sty.n_workers = 2
    sty.file_save_mode = :append
    sty.save_folder = mktempdir()
    sty.n_simulations = 2
    sty.release_workers_after_run = false


    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    act.ActorPedestal.density_match = :ne_line
    ini.core_profiles.ne_setting = :greenwald_fraction
    ini.core_profiles.ne_value = 0.2 ↔ [0.2, 1.0]

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

    @info "Running study... "
    FUSE.run(study)

    @testset "inis and acts" begin
        ini_list = [rand(ini) for _ in 1:2]
        act_list = [rand(act) for _ in 1:2]

        study = FUSE.StudyDatabaseGenerator(sty, ini_list, act_list)
        study.workflow = workflow_DatabaseGenerator

        @info "Running study for predefined `inis` and `acts`..."
        FUSE.run(study)
    end

    # Manual release of workers
    workers_list = Distributed.workers()
    if !isempty(workers_list)
        Distributed.rmprocs(workers_list)
        @info "released $(length(workers_list)) workers: $workers_list"
    else
        @info "no workers to release"
    end
end