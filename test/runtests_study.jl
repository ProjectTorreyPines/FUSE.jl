using FUSE
using Test
using Distributed

@testset
    sty,act = FUSE.study_parameters(:DatabaseGenerator);
    sty.server = "localhost"
    sty.n_workers = 2

    sty.file_save_mode = :append
    sty.save_folder = mkdtempdir()
    sty.n_simulations = 2
    ini, act = FUSE.case_parameters(:ITER;init_from=:scalars)
    ini.core_profiles.ne_setting = :greenwald_fraction
    ini.core_profiles.ne_value = 0.2 ↔ [0.2,1.0]
    study = FUSE.StudyDatabaseGenerator(sty,ini, act); # it is possible to pass in keyword arguments to sty

    @everywhere import FUSE
    @everywhere function workflow_DatabaseGenerator(dd::FUSE.IMAS.dd, ini::FUSE.ParametersAllInits, act::FUSE.ParametersAllActors)
        FUSE.init(dd, ini, act)
        return nothing
    end
    FUSE.run(study)
end