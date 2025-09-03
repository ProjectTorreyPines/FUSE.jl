using Revise
using FUSE
using IMAS;
using Plots;
FUSE.logging(Logging.Info; actors=Logging.Info);
# sty is the act equivalent for a study, it has common parameters like server and n_workers but also study dependent parameters like TGLF saturation rules
sty,act = FUSE.study_parameters(:TGLFdb);
sty
# Interacting with sty
sty.server = "saga" #localhost #saga
sty.n_workers = 24

sty.database_folder = "/mnt/beegfs/users/neisert/ODSs/d3d_negDcake_oak"

#mastu_intersect_steady "#d3d_negDcake_oak"#iri/166066"#d3d" #201_all_nstx #test_nstx

sty.save_folder = "/mnt/beegfs/users/neisert/ODSs/d3d_negDcake_oak/outputs_edgerun_rho_0p1_0p021_0p99_rot"

sty.sat_rules = missing #[:sat1,:sat2,:sat3] #study specific parameters
# It's also possible to run with a custom tglfnn model, set sty.custom_tglf_models
sty.custom_tglf_models =  ["sat1_em_d3d_negD","sat1_em_d3d_azf-1","sat3_em_d3d_azf+1", "sat3_em_d3d_azf-1"]

#"sat0quench_em_d3d_azf+1", "sat0_em_d3d", "sat1_em_d3d_negD", "sat2_em_d3d_negD", "sat3_em_d3d_negD", "sat1_em_d3d_azf+1", "sat1_em_d3d_azf-1", "sat2_em_d3d_azf+1", "sat2_em_d3d_azf-1",  "sat3_em_d3d_azf+1", "sat3_em_d3d_azf-1"]

#["sat2_em_nstx_azf-1"]
#,"sat3_em_nstx_azf-1"]
#["sat0quench_em_nstx_azf+1","sat1_em_nstx_azf-1","sat1geo_em_nstx_azf-1"]#, "sat2_em_nstx_azf-1","sat3_em_nstx_azf-1"]


#["sat0quench_es_d3d_azf+1", "sat0quench_em_d3d_azf+1", "sat0_es_d3d", "sat0_em_d3d", "sat1_es_d3d_azf+1", "sat1_es_d3d_azf-1", "sat1_em_d3d_azf+1", "sat1_em_d3d_azf-1", "sat2_es_d3d_azf+1", "sat2_es_d3d_azf-1", "sat2_em_d3d_azf+1", "sat2_em_d3d_azf-1", "sat3_es_d3d_azf+1", "sat3_es_d3d_azf-1", "sat3_em_d3d_azf+1", "sat3_em_d3d_azf-1"]

#["sat1_em_d3d_azf+1", "sat2_em_d3d_azf-1", "sat3_em_d3d_azf-1", "sat3_em_d3d_azf+1", "sat2_em_d3d_azf+1"]

#["sat1_em_mastu_azf-1"] #,"sat1_es_d3d_azf-1"]  #, "sat3_em_d3d_azf-1"] #, "sat1_em_d3d_azf+1", "sat3_em_d3d_azf+1"]

# All DIII-D
# ["sat0quench_es_d3d_azf+1", "sat0quench_em_d3d_azf+1", "sat0_es_d3d", "sat0_em_d3d", "sat1_es_d3d_azf+1", "sat1_es_d3d_azf-1", "sat1_em_d3d_azf+1", "sat1_em_d3d_azf-1", "sat2_es_d3d_azf+1", "sat2_es_d3d_azf-1", "sat2_em_d3d_azf+1", "sat2_em_d3d_azf-1", "sat3_es_d3d_azf+1", "sat3_es_d3d_azf-1", "sat3_em_d3d_azf+1", "sat3_em_d3d_azf-1"]

# All MAST-U
#["sat1_es_d3d+mastu_azf+1", "sat0quench_em_d3d+mastu_azf+1","sat0quench_es_mastu_azf+1","sat0quench_em_mastu_azf+1","sat0_es_mastu_azf+1", "sat0_em_mastu_azf+1", "sat1_es_mastu_azf+1", "sat1_em_mastu_azf+1","sat1_es_mastu_azf-1", "sat1_em_mastu_azf-1", "sat2_es_mastu_azf+1", "sat2_es_mastu_azf-1", "sat2_em_mastu_azf+1", "sat2_em_mastu_azf-1", "sat3_es_mastu_azf+1", "sat3_es_mastu_azf-1", "sat3_em_mastu_azf+1", "sat3_em_mastu_azf-1"]

# All negD
#["sat0_em_d3d_negD", "sat1_em_d3d_negD", "sat2_em_d3d_negD", "sat3_em_d3d_negD"]

# FPP
#["sat1_em_fpp", "sat1_em_fpp_iter", "sat1_em_fpp_d3d"] 

sty.file_save_mode = :overwrite #overwrite or #safe_write
sty.release_workers_after_run = true # this is the default behavior and releases workers after running the study
sty
study = FUSE.StudyTGLFdb(sty, act); # it is possible to pass in keyword arguments to sty
using Distributed
@everywhere import FUSE
#change density evolution in TGLFdb.jl
study.act.ActorFluxMatcher.evolve_rotation = :flux_match #fixed nstx
study.act.ActorFluxMatcher.rho_transport = sort(vcat(0.1:0.02:0.90, 0.91:0.01:0.99)) #watch step and init
study.act.ActorFluxMatcher
FUSE.run(study);
FUSE.analyze(study);