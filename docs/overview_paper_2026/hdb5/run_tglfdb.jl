using FUSE

using FUSE

sty = FUSE.study_parameters(:TGLFdb)   # or the exact symbol your install exposes
sty.server = "localhost"
sty.n_workers = 127
sty.release_workers_after_run = true
sty.save_folder = "tglfdb_results"

study = FUSE.StudyTGLFdb(sty)
FUSE.run(study)
