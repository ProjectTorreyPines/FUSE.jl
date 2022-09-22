# SPARC and ARC use cases


```@julia
using Revise
using FUSE
using Plots; gr();
FUSE.logging(Logging.Info);
```

### Initialize `dd`, `ini`, `act` from SPARC use case


```@julia
dd, ini, act = FUSE.init(:SPARC, do_plot=true);
```

### Plot FUSE build on top of SPARC drawing



```@julia
CAD = FUSE.TraceCAD(:SPARC)
plot(CAD, size=(900,900))
plot!(dd.equilibrium, cx=true, color=:red)
plot!(dd.build, wireframe=true, linewidth=2, color=:black)
```

--------------

### Initialize `dd`, `ini`, `act` from ARC use case


```@julia
dd, ini, act = FUSE.init(:ARC, do_plot=true);
```

### Plot FUSE build on top of ARC drawinga



```@julia
CAD = FUSE.TraceCAD(:ARC)
plot(CAD, size=(900,900))
plot!(dd.equilibrium, cx=true, color=:red)
plot!(dd.build, wireframe=true, linewidth=2, color=:black)
```
