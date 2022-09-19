# Blanket


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Initialize the ITER case case
[ITER case documentation](https://fuse.help/cases.html#ITER)


```@julia
dd, ini, act = FUSE.init(:ITER, init_from=:ods, do_plot=true);
```

### Run Actors that will be needed for the blanket actor


```@julia
FUSE.ActorEquilibriumTransport(dd, act)
FUSE.ActorCXbuild(dd, act)
FUSE.ActorNeutronics(dd, act; do_plot=true);
```

### Running the simple blanket actor
[ActorBlanket documentation](https://fuse.help/actors.html#Blanket)


```@julia
dd.build.structure
FUSE.ActorBlanket(dd, act);
dd.blanket
```

### Running the blanket actor
[ActorBlanket documentation](https://fuse.help/actors.html#Blanket)


```@julia
FUSE.ActorBlanket(dd, act)
dd.blanket
```


```@julia
act
```
