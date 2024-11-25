using FUSE
using Test
using Distributed
using ProgressMeter

@testset "study" begin
    sty, act = FUSE.study_parameters(:DatabaseGenerator)
    sty.server = "localhost"
    sty.n_workers = 2
    sty.file_save_mode = :append
    sty.save_folder = mktempdir()
    sty.n_simulations = 2

    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    act.ActorPedestal.density_match=:ne_line
    ini.core_profiles.ne_setting = :greenwald_fraction
    ini.core_profiles.ne_value = 0.2 ↔ [0.2, 1.0]

    @everywhere import FUSE

    @everywhere function workflow_DatabaseGenerator(dd::FUSE.IMAS.dd, ini::FUSE.ParametersAllInits, act::FUSE.ParametersAllActors)
        FUSE.init(dd, ini, act)
        return nothing
    end

    study = FUSE.StudyDatabaseGenerator(sty, ini, act)
    study.workflow = workflow_DatabaseGenerator

    FUSE.run(study)
end