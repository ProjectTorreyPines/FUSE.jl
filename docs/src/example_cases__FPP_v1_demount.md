# GA's FPP use case


```@julia
using Revise
using FUSE
using Plots; gr();
FUSE.logging(Logging.Info);
```

### Get `ini` and `act` for FPP use case starting from scalars and taking the result from STEP in OMFIT (:ods)
* [ini documentation](https://fuse.help/ini.html)
* [act documentation](https://fuse.help/act.html)


```@julia
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars);
ini_ods, act_ods = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:ods);
```

### Customize `ini` (or `act`) parameters


```@julia
# See the content of ini
# ini.equilibrium.ζ = 0.15 # bump up squareness (requires CHEASE)
ini
```


```@julia
# See the content of act
# act.ActorEquilibrium.model = :CHEASE # use CHEASE?
act
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
# look at what parameters for example ActorStationaryPlasma
display(act.ActorEquilibriumTransport)

# look at the details some of these parameters
display(act.ActorStationaryPlasma)


```


```@julia
# run the actor
FUSE.ActorStationaryPlasma(dd, act; do_plot=true);
```

### Initialize FPP case from STEP


```@julia
dd_ods = IMAS.dd()
FUSE.init(dd_ods, ini_ods, act_ods; do_plot=true);
```


```@julia
dd_ods = IMAS.dd()
FUSE.init(dd_ods, ini_ods, act_ods; do_plot=true);
dd_ods.core_profiles.global_quantities.ejima = [0.4]

# Running whole facility without transport as this was carefuly done in STEP to obtain the ITB scenario

# Solvev into CHEASE to form the x-point from an equilibrium without x-points
act_ods.ActorEquilibrium.model = :Solovev
FUSE.ActorEquilibrium(dd_ods,act_ods)
act_ods.ActorEquilibrium.model = :CHEASE
FUSE.ActorEquilibrium(dd_ods,act_ods)

FUSE.ActorHFSsizing(dd_ods,act_ods)
FUSE.ActorLFSsizing(dd_ods,act_ods)
FUSE.ActorCXbuild(dd_ods, act_ods)
FUSE.ActorNeutronics(dd_ods,act_ods)
FUSE.ActorBlanket(dd_ods,act_ods)
FUSE.ActorDivertors(dd_ods,act_ods)
FUSE.ActorBalanceOfPlant(dd_ods,act_ods)
plot(dd_ods.build)
display(plot!(dd_ods.equilibrium,cx=true))
plot(dd_ods.core_profiles)
println("Total net electricity out of this ITB FPP case from running STEP transport-equilibrium $(round(dd_ods.balance_of_plant.power_electric_net[1]/1e6)) MWe")
```
