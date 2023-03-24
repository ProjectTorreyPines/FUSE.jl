# Balance of plant


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info; actors=Logging.Debug);
```

### Initialize FPP v1_demount case
[FPP v1 demount case documentation](https://fuse.help/cases.html#FPP)


```@julia
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars);
dd = FUSE.init(ini, act; do_plot=false);
```

### Run Actors that will be needed for balance of plant


```@julia
FUSE.ActorNeutronics(dd, act; do_plot=false)
FUSE.ActorDivertors(dd, act)
FUSE.ActorBlanket(dd, act);
```

### Running the simple brayton cycle
Run the balance of plant actor, with for the simple case of "brayton only", this is a generalized power cycle which does not optimize itself depending on the configuration


```@julia
empty!(dd.balance_of_plant)
act.ActorThermalCycle.power_cycle_type=:brayton_only
FUSE.ActorBalanceOfPlant(dd, act; do_plot=true)
display(IMAS.freeze(dd.balance_of_plant))
```

### Running the complex_brayton file
[ActorBalanceOfPlant documentation](https://fuse.help/actors.html#BalanceOfPlant)
Run the balance of plant with the model "complex_brayton". This configuration modifies the component order and operating temperatures to minimize the heat waste. For this case, the complex_brayton model has a thermal effeciency of 45%. 
This model relies on high operating temperatures


```@julia
empty!(dd.balance_of_plant)
act.ActorThermalCycle.power_cycle_type=:complex_brayton
FUSE.ActorBalanceOfPlant(dd, act; do_plot = true);
display(IMAS.freeze(dd.balance_of_plant))
```

## Power scan for different power cycles

Basic fusion power scan that only runs the ActorBalanceOfPlant, to see what is the break-even point for the different power cycles that we support. This scan knows nothing about the plasma, it only cares about the thermal power that is extracted at the blanket and divertors.


```@julia
act = FUSE.ParametersActors()
dd = IMAS.dd()

# minimal initialization to run ActorBalanceOfPlant
dd.core_profiles.time = [0.0]
resize!(dd.blanket.module, 1)
resize!(dd.blanket.module[1].time_slice)
dd.blanket.module[1].time_slice[].power_thermal_extracted = 0.0
resize!(dd.divertors.divertor, 1)
@ddtime dd.divertors.divertor[1].power_incident.data = 0.0

# assume 80% on blanket and rest on divertors
blanket_div_ratio = 0.8

# extent of the scan
POW = LinRange(6, 9, 100)

# run and plot
p=plot()
for cycle_type in [:brayton_only, :complex_brayton]#, :rankine_only
    
    Penet = Float64[]
    Ppump = Float64[]
    act.ActorThermalCycle.power_cycle_type=cycle_type
    for pow in POW
        dd.blanket.module[1].time_slice[].power_thermal_extracted = 10^(pow) * blanket_div_ratio
        @ddtime dd.divertors.divertor[1].power_incident.data = 10^(pow) * (1.0 - blanket_div_ratio)
        FUSE.ActorBalanceOfPlant(dd, act)
        push!(Penet, dd.balance_of_plant.power_electric_net[1])
        push!(Ppump, findfirst(:pumping, dd.balance_of_plant.power_electric_plant_operation.system).power[1])
    end

    pow0 = POW[argmin(abs.(Penet))]
    plot!(p, 10.0 .^ POW, Penet, xscale=:log10, label="$cycle_type $(round(10^pow0/1E6)) MW",
        xlabel="Power thermal extracted [W]", ylabel="Power_electric net [W]",legend=:topleft)
    #plot!(p, 10.0 .^ POW, -Ppump, primary=false)
    vline!(p, [10^pow0], label="", primary=false, ls=:dash)
end
hline!(p, [0.0], label="", color=:black)
```
