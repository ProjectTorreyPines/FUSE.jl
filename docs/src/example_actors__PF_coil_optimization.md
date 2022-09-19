# Poloidal field coil optimization (ITER)


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Initialization of `dd`, `ini`, `act` for the ITER use case


```@julia
dd, ini, act = FUSE.init(:ITER, init_from=:ods, do_plot=true);
```

### Run the SteadyStateCurrent actor in order to estimate how much ohmic current will be required during flattop
[SteadyStateCurrent actor](https://fuse.help/actors.html#SteadyStateCurrent)



```@julia
FUSE.ActorSteadyStateCurrent(dd, act);
```

### Estimate how much flux is required during start-up
[FluxSwing actor](https://fuse.help/actors.html#FluxSwing)




```@julia
FUSE.ActorFluxSwing(dd, act; operate_at_j_crit=true, j_tolerance=0.0)
# critial_j shoudl be the same as max_j
dd.build.oh
# Estimated flattop durration at maximum oh current operation suggests that ITER can run for 3500 seconds! 
```

### Run PF coils optimization and plot


```@julia
FUSE.init_pf_active(dd, ini, act)
#plot(dd.pf_active)
#act.ActorPFcoilsOpt[:optimization_scheme]
FUSE.ActorPFcoilsOpt(dd, act; optimization_scheme=:currents, do_plot=true);
```
