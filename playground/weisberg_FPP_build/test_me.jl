using Revise
using FUSE
using IMAS
using Plots;
FUSE.logging(Logging.Info; actors=Logging.Debug);

dd = IMAS.dd();
ini, act = FUSE.case_parameters(:FPP, version=:v2);

act.ActorHFSsizing.error_on_technology = false
act.ActorHFSsizing.error_on_performance = false

FUSE.init(dd, ini, act; do_plot=false);

#@profview FUSE.ActorSteadyStateCurrent(dd, act);
@profview FUSE.ActorHFSsizing(dd, act);
