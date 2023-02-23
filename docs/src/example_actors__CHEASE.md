# CHEASE


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Initialize the ITER case from ODS
[ITER case documentation](https://fuse.help/cases.html#ITER)


```@julia
dd, ini, act = FUSE.init(:ITER, init_from=:ods);
```

#### Let's run the ActorCHEASE and compare with original equlibrium

[CHEASE actor](https://fuse.help/actors.html#CHEASE) documentation

NOTE: CHEASE is a fixed boundary equilibrium solver. Extension of the magnetic field calculation in the vacuum region (including X-points) can be triggered by setting `act.ActorCHEASE.free_boundary = true`. This is the default behavior in FUSE, since diverted plasma information is needed to proceed with build and divertor actors.


```@julia
dd, ini, act = FUSE.init(:ITER, init_from=:ods);
eq_plot = plot(dd.equilibrium, label="original eq")
act.ActorCHEASE.free_boundary = true
actor = FUSE.ActorCHEASE(dd, act);
display(plot!(eq_plot, dd.equilibrium, label="CHEASE w/ vacuum fields"))
```


```@julia
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:ods);
dd = IMAS.dd()
FUSE.init(dd, ini, act; do_plot=true);
```

### Runnning CHEASE on all configurations of the signs of Ip and Bt


```@julia
sample_dir = abspath(joinpath(@__DIR__, "../../sample"))
files = [joinpath(sample_dir, file) for file in readdir(sample_dir) if startswith(file, "g") && endswith(file, ".json")]
act = FUSE.ParametersActors()

for file in files
    dd = IMAS.json2imas(file)
    plot(dd.equilibrium, label="$(split(split(file,"/")[end],".")[1]) EQ")
    FUSE.ActorCHEASE(dd, act)
    label = "CHEASE EQ Bt = $(round(@ddtime(dd.equilibrium.vacuum_toroidal_field.b0),digits=2)) [T], Ip = $(round(dd.equilibrium.time_slice[].global_quantities.ip/1e6,digits=2)) [MA]"
    display(plot!(dd.equilibrium, label=label))
end
```
