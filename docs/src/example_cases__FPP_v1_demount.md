# GA's FPP use case


```@julia
using Revise
using FUSE
using Plots; gr();
FUSE.logging(Logging.Info);
```

### Get `ini` and `act` for FPP use case
* [ini documentation](https://fuse.help/ini.html)
* [act documentation](https://fuse.help/act.html)


```@julia
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars);
```

### Customize `ini` (or `act`) parameters
* Changing Zeff from 1.1 to 2.0 will improve confinement significantly due to the pedestal increase!


```@julia
ini.core_profiles.zeff = 2.0
ini.core_profiles # see the content of ini.core_profiles
```

### Initialize `dd` based on the `ini` and `act` parameters

* [Data structure (dd) documentation](https://fuse.help/dd.html)


```@julia
dd = IMAS.dd()
FUSE.init(dd, ini, act; do_plot=true);
```

### Run the coupled equilibrium-transport actor 

* [Equilibrium Transport actor documentation](https://fuse.help/actors.html#EquilibriumTransport)


```@julia
# look at what parameters for example ActorEquilibriumTransport and ActorTauenn use
display(act.ActorEquilibriumTransport)
display(act.ActorTauenn)

# look at the details some of these parameters
display(act.ActorTauenn[:temp_pedestal_ratio])

# modify parameters
act.ActorTauenn.temp_pedestal_ratio = 0.9;
```


```@julia
# run the actor
FUSE.ActorEquilibriumTransport(dd, act; do_plot=true);
```
