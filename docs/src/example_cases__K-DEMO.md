# K-DEMO use case


```@julia
using FUSE
using Plots;
FUSE.logging(Logging.Info; actors=Logging.Info);
```


```@julia
ini, act = FUSE.case_parameters(:KDEMO)
display(ini.equilibrium)
```


```@julia
@show ini.oh.technology, ini.tf.technology, ini.pf_active.technology
# a note that this is ITER-like LTS
```


```@julia
dd_LTS = FUSE.init(ini, act);
```


```@julia
dd = deepcopy(dd_LTS)

act.ActorStabilityLimits.raise_on_breach = false
act.ActorHFSsizing.error_on_technology = false

FUSE.ActorWholeFacility(dd, act);
```


```@julia
# Make magnets HTS (ideally we get K-DEMO LTS magnet techonolgy parameters)
ini.oh.technology = :HTS
ini.tf.technology = :HTS
ini.pf_active.technology = :HTS

dd_HTS = FUSE.init(ini, act)
act.ActorPFcoilsOpt.do_plot = true

act.ActorHFSsizing.verbose = true

FUSE.ActorWholeFacility(dd_HTS, act);
```


```@julia
FUSE.digest(dd_HTS)
```


```@julia
# Core transport according to turbulent and neoclassical fluxes with EPED H-mode at the edge

act.ActorFluxMatcher.max_iterations = 300
act.ActorFluxMatcher.evolve_pedestal = true
act.ActorFluxMatcher.verbose = true
act.ActorTGLF.user_specified_model = "sat1_em_iter"
FUSE.ActorFluxMatcher(dd_HTS, act);
```


```@julia
FUSE.digest(dd_HTS, "K DEMO using FUSE WholeFacility")
```
