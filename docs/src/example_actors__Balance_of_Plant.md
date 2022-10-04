# Balance of plant


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Initialize FPP v1_demount case
[FPP v1 demount case documentation](https://fuse.help/cases.html#FPP)


```@julia
dd, ini, act = FUSE.init(:FPP, version=:v1_demount, init_from=:ods, do_plot=false);
```

### Run Actors that will be needed for balance of plant


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

### Running the divertors actor
[ActorDivertors documentation](https://fuse.help/actors.html#Divertors)


```@julia
FUSE.ActorDivertors(dd, act)
dd.divertors
```

### Running the balance of plant actor
[ActorBalanceOfPlant documentation](https://fuse.help/actors.html#BalanceOfPlant)


```@julia
FUSE.ActorBalanceOfPlant(dd, act);

println("The net electrical power to the grid is $(round(dd.balance_of_plant.power_electric_net[end]/1e6,digits=1)) [MWe] \nWith Qplant = $(round(dd.balance_of_plant.Q_plant[end],digits=2)) \n")
display(dd.balance_of_plant)

```
