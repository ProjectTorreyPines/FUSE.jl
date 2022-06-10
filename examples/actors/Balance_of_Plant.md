# Balance of plant


```julia
using Revise
using FUSE
using Plots; gr();
global_logger(FUSE.logger);
```

### Initialize FPP v1_demount case
[FPP v1 demount case documentation](https://fuse.help/cases.html#FPP)


```julia
dd, ini, act = FUSE.init(:FPP, version=:v1_demount, init_from=:ods, do_plot=false);
```

    WARNING: both ImageMetadata and ImageAxes export "data"; uses of it in module Images must be qualified
    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mdd.dataset_description was skipped in IMAS data dictionary
    [33m[1m└ [22m[39m[90m@ IMASDD ~/.julia/dev/IMASDD/src/data.jl:1067[39m


### Run Actors that will be needed for balance of plant


```julia
FUSE.ActorEquilibriumTransport(dd, act)
FUSE.ActorHFSsizing(dd, act)
FUSE.ActorLFSsizing(dd, act)
FUSE.ActorCXbuild(dd, act)
FUSE.ActorNeutronics(dd, act; do_plot=true);
```


    
![svg](Balance_of_Plant_files/Balance_of_Plant_5_0.svg)
    


### Running the simple blanket actor
[ActorBlanket documentation](https://fuse.help/actors.html#Blanket)


```julia
dd.build.structure
FUSE.ActorBlanket(dd,act);
dd.blanket
```




    [0m[1mblanket[22m
    ├─ [0m[1mmodule[22m
    │  ├─ [0m[1m1[22m
    │  │  ├─ [0mname[31m ➡ [39m[35m"LFS blanket"[39m
    │  │  └─ [0m[1mtime_slice[22m
    │  │     └─ [0m[1m1[22m
    │  │        ├─ [0mpower_incident_neutrons[31m ➡ [39m[31m6.11188e+08[39m
    │  │        ├─ [0mpower_incident_radiated[31m ➡ [39m[31m0[39m
    │  │        ├─ [0mpower_thermal_extracted[31m ➡ [39m[31m7.33425e+08[39m
    │  │        ├─ [0mpower_thermal_neutrons[31m ➡ [39m[31m7.33425e+08[39m
    │  │        ├─ [0mpower_thermal_radiated[31m ➡ [39m[31m0[39m
    │  │        └─ [0mtritium_breeding_ratio[31m ➡ [39m[31m1[39m
    │  └─ [0m[1m2[22m
    │     ├─ [0mname[31m ➡ [39m[35m"HFS blanket"[39m
    │     └─ [0m[1mtime_slice[22m
    │        └─ [0m[1m1[22m
    │           ├─ [0mpower_incident_neutrons[31m ➡ [39m[31m1.69757e+08[39m
    │           ├─ [0mpower_incident_radiated[31m ➡ [39m[31m0[39m
    │           ├─ [0mpower_thermal_extracted[31m ➡ [39m[31m2.03708e+08[39m
    │           ├─ [0mpower_thermal_neutrons[31m ➡ [39m[31m2.03708e+08[39m
    │           ├─ [0mpower_thermal_radiated[31m ➡ [39m[31m0[39m
    │           └─ [0mtritium_breeding_ratio[31m ➡ [39m[31m1[39m
    ├─ [0mtime[31m ➡ [39m[32m[1][39m
    └─ [0mtritium_breeding_ratio[31m ➡ [39m[32m[0.928583][39m




### Running the divertors actor
[ActorDivertors documentation](https://fuse.help/actors.html#Divertors)


```julia
FUSE.ActorDivertors(dd,act)
dd.divertors
```




    [0m[1mdivertors[22m
    ├─ [0m[1mdivertor[22m
    │  ├─ [0m[1m1[22m
    │  │  ├─ [0mname[31m ➡ [39m[35m"Upper divertor"[39m
    │  │  ├─ [0m[1mpower_incident[22m
    │  │  │  ├─ [0mdata[31m ➡ [39m[32m[1.47564e+08][39m
    │  │  │  └─ [0mtime[31m ➡ [39m[32m[1][39m
    │  │  ├─ [0m[1mpower_neutrals[22m
    │  │  │  ├─ [0mdata[31m ➡ [39m[32m[0][39m
    │  │  │  └─ [0mtime[31m ➡ [39m[32m[1][39m
    │  │  ├─ [0m[1mpower_radiated[22m
    │  │  │  ├─ [0mdata[31m ➡ [39m[32m[0][39m
    │  │  │  └─ [0mtime[31m ➡ [39m[32m[1][39m
    │  │  ├─ [0m[1mpower_recombination_neutrals[22m
    │  │  │  ├─ [0mdata[31m ➡ [39m[32m[0][39m
    │  │  │  └─ [0mtime[31m ➡ [39m[32m[1][39m
    │  │  └─ [0m[1mpower_thermal_extracted[22m
    │  │     ├─ [0mdata[31m ➡ [39m[32m[1.47564e+08][39m
    │  │     └─ [0mtime[31m ➡ [39m[32m[1][39m
    │  └─ [0m[1m2[22m
    │     ├─ [0mname[31m ➡ [39m[35m"Lower divertor"[39m
    │     ├─ [0m[1mpower_incident[22m
    │     │  ├─ [0mdata[31m ➡ [39m[32m[1.47564e+08][39m
    │     │  └─ [0mtime[31m ➡ [39m[32m[1][39m
    │     ├─ [0m[1mpower_neutrals[22m
    │     │  ├─ [0mdata[31m ➡ [39m[32m[0][39m
    │     │  └─ [0mtime[31m ➡ [39m[32m[1][39m
    │     ├─ [0m[1mpower_radiated[22m
    │     │  ├─ [0mdata[31m ➡ [39m[32m[0][39m
    │     │  └─ [0mtime[31m ➡ [39m[32m[1][39m
    │     ├─ [0m[1mpower_recombination_neutrals[22m
    │     │  ├─ [0mdata[31m ➡ [39m[32m[0][39m
    │     │  └─ [0mtime[31m ➡ [39m[32m[1][39m
    │     └─ [0m[1mpower_thermal_extracted[22m
    │        ├─ [0mdata[31m ➡ [39m[32m[1.47564e+08][39m
    │        └─ [0mtime[31m ➡ [39m[32m[1][39m
    └─ [0mtime[31m ➡ [39m[32m[1][39m




### Running the balance of plant actor
[ActorBalanceOfPlant documentation](https://fuse.help/actors.html#BalanceOfPlant)


```julia
FUSE.ActorBalanceOfPlant(dd,act);

println("The net electrical power to the grid is $(round(dd.balance_of_plant.power_electric_net[end]/1e6,digits=1)) [MWe] \nWith Qplant = $(round(dd.balance_of_plant.Q_plant[end],digits=2)) \n")
display(dd.balance_of_plant)

```

    The net electrical power to the grid is 242.9 [MWe] 
    With Qplant = 1.97 
    



    [0m[1mbalance_of_plant[22m
    ├─ [0mQ_plant[31m ➡ [39m[34mFunction[39m
    ├─ [0mpower_electric_net[31m ➡ [39m[34mFunction[39m
    ├─ [0m[1mpower_electric_plant_operation[22m
    │  ├─ [0m[1msystem[22m
    │  │  ├─ [0m[1m1[22m
    │  │  │  ├─ [0mindex[31m ➡ [39m[33m1[39m
    │  │  │  ├─ [0mname[31m ➡ [39m[35m"H&CD"[39m
    │  │  │  ├─ [0mpower[31m ➡ [39m[32m100-element Vector{Float64}[39m
    │  │  │  └─ [0m[1msubsystem[22m
    │  │  │     └─ [0m[1m1[22m
    │  │  │        ├─ [0mindex[31m ➡ [39m[33m1[39m
    │  │  │        ├─ [0mname[31m ➡ [39m[35m"ec_launchers"[39m
    │  │  │        └─ [0mpower[31m ➡ [39m[32m100-element Vector{Float64}[39m
    │  │  ├─ [0m[1m2[22m
    │  │  │  ├─ [0mindex[31m ➡ [39m[33m2[39m
    │  │  │  ├─ [0mname[31m ➡ [39m[35m"cryostat"[39m
    │  │  │  └─ [0mpower[31m ➡ [39m[32m100-element Vector{Float64}[39m
    │  │  ├─ [0m[1m3[22m
    │  │  │  ├─ [0mindex[31m ➡ [39m[33m3[39m
    │  │  │  ├─ [0mname[31m ➡ [39m[35m"tritium_handling"[39m
    │  │  │  └─ [0mpower[31m ➡ [39m[32m100-element Vector{Float64}[39m
    │  │  ├─ [0m[1m4[22m
    │  │  │  ├─ [0mindex[31m ➡ [39m[33m4[39m
    │  │  │  ├─ [0mname[31m ➡ [39m[35m"pumping"[39m
    │  │  │  └─ [0mpower[31m ➡ [39m[32m100-element Vector{Float64}[39m
    │  │  └─ [0m[1m5[22m
    │  │     ├─ [0mindex[31m ➡ [39m[33m5[39m
    │  │     ├─ [0mname[31m ➡ [39m[35m"pf_active"[39m
    │  │     └─ [0mpower[31m ➡ [39m[32m100-element Vector{Float64}[39m
    │  └─ [0mtotal_power[31m ➡ [39m[34mFunction[39m
    ├─ [0m[1mthermal_cycle[22m
    │  ├─ [0mpower_electric_generated[31m ➡ [39m[34mFunction[39m
    │  ├─ [0m[1msystem[22m
    │  │  ├─ [0m[1m1[22m
    │  │  │  ├─ [0mindex[31m ➡ [39m[33m1[39m
    │  │  │  ├─ [0mname[31m ➡ [39m[35m"blanket"[39m
    │  │  │  └─ [0mpower_in[31m ➡ [39m[32m100-element Vector{Float64}[39m
    │  │  └─ [0m[1m2[22m
    │  │     ├─ [0mindex[31m ➡ [39m[33m2[39m
    │  │     ├─ [0mname[31m ➡ [39m[35m"divertors"[39m
    │  │     └─ [0mpower_in[31m ➡ [39m[32m100-element Vector{Float64}[39m
    │  └─ [0mthermal_electric_conversion_efficiency[31m ➡ [39m[32m100-element Vector{Float64}[39m
    └─ [0mtime[31m ➡ [39m[32m100-element Vector{Float64}[39m




```julia

```
