# Neutronics


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Initialize FPP v1_demount case
[FPP v1 demount case documentation](https://fuse.help/cases.html#FPP)


```@julia
dd, ini, act = FUSE.init(:FPP, version=:v1_demount, init_from=:scalars, do_plot=false);
```

### Run the Neutronics actor


```@julia
FUSE.ActorNeutronics(dd, act; do_plot=true);
```

### Plot neutron wall loading


```@julia
plot(dd.neutronics.time_slice[].wall_loading; cx=true)
```


```@julia
plot(dd.neutronics.time_slice[].wall_loading; cx=false)
```
