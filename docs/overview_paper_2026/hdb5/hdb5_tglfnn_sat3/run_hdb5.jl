using FUSE

sty = FUSE.study_parameters(:HDB5Validation)

sty.server = "localhost"

# Use most or all physical cores as workers, but remember the master is also a process.
# 127 is the maximal "fill the node" choice on a 128-core CPU node.
sty.n_workers = 127

sty.release_workers_after_run = true
sty.save_folder = "hdb5_tglfnn_sat3"
sty.tokamak = :all
sty.n_samples_per_tokamak = 0  # all cases; set smaller for testing

_, act = FUSE.case_parameters(:HDB5; tokamak=:JET)
#act.ActorTGLF.model=:TJLF # for full TJLF...  painfully slow
act.ActorTGLF.tglfnn_model = "sat3_em_d3d_azf-1_withnegD"

study = FUSE.StudyHDB5Validation(sty, act)
FUSE.run(study)
