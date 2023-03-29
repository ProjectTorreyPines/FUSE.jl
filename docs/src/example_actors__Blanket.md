# Blanket


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Initialize case


```@julia
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars);
dd = FUSE.init(ini, act; do_plot=false);
```

### Run Neutronics actors


```@julia
FUSE.ActorNeutronics(dd, act; do_plot=true);
```

### Running the blanket actor
[ActorBlanket documentation](https://fuse.help/actors.html#Blanket)


```@julia
dd.build.structure
act.ActorBlanket.minimum_first_wall_thickness = 0.02
FUSE.ActorBlanket(dd, act, verbose=true);
dd.blanket
```
