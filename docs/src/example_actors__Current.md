# Current evolution


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info; actors=Logging.Debug);
```


```@julia
# Load TRANSP data at 2.91 s
import QED
file_0 = joinpath(dirname(pathof(QED)), "..","sample", "ods_163303Z29-2910.json")
dd = IMAS.json2imas(file_0; verbose=false)

act=FUSE.ParametersActors()

# initialize actor
actor = FUSE.ActorQED(dd, act)

# evolve current
#for k in 1:3
#    FUSE.step(actor, 0.1, 100, resume=true)
#    FUSE.finalize(actor)
#    dd.global_time = dd.equilibrium.time[end]
#end
```

    [ Info: QED



    MethodError: no method matching _step(::FUSE.ActorQED{Float64, Float64})
    
    Closest candidates are:
      _step(::FUSE.ActorQED, ::Any, ::Any)
       @ FUSE ~/.julia/dev/FUSE/src/actors/current/qed_actor.jl:46
      _step(::FUSE.ActorQED, ::Any, ::Any, ::Any)
       @ FUSE ~/.julia/dev/FUSE/src/actors/current/qed_actor.jl:46
      _step(::FUSE.ActorQED, ::Any, ::Any, ::Any, ::Any; resume)
       @ FUSE ~/.julia/dev/FUSE/src/actors/current/qed_actor.jl:46
      ...


    

    Stacktrace:

     [1] macro expansion

       @ ~/.julia/dev/FUSE/src/actors/abstract_actors.jl:41 [inlined]

     [2] macro expansion

       @ ~/.julia/packages/TimerOutputs/RsWnF/src/TimerOutput.jl:237 [inlined]

     [3] step(::FUSE.ActorQED{Float64, Float64}; kw::Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}})

       @ FUSE ~/.julia/dev/FUSE/src/actors/abstract_actors.jl:38

     [4] step(::FUSE.ActorQED{Float64, Float64})

       @ FUSE ~/.julia/dev/FUSE/src/actors/abstract_actors.jl:33

     [5] FUSE.ActorQED(dd::IMASdd.dd{Float64}, act::FUSE.ParametersActors{Float64}; kw::Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}})

       @ FUSE ~/.julia/dev/FUSE/src/actors/current/qed_actor.jl:31

     [6] FUSE.ActorQED(dd::IMASdd.dd{Float64}, act::FUSE.ParametersActors{Float64})

       @ FUSE ~/.julia/dev/FUSE/src/actors/current/qed_actor.jl:29

     [7] top-level scope

       @ In[27]:9



```@julia
dd
```




    dd
    ├─ balance_of_plant
    │  ├─ Q_plant ➡ Function
    │  ├─ power_electric_net ➡ Function [W]
    │  ├─ power_electric_plant_operation
    │  │  └─ total_power ➡ Function [W]
    │  └─ thermal_cycle
    │     ├─ power_electric_generated ➡ Function [W]
    │     └─ total_useful_heat_power ➡ Function [W]
    ├─ build
    │  └─ tf
    │     ├─ ripple ➡ Function
    │     └─ wedge_thickness ➡ Function [m]
    ├─ core_profiles
    │  ├─ profiles_1d
    │  │  └─ 1
    │  │     ├─ conductivity_parallel ➡ 101-element Vector{Float64} [ohm^-1.m^-1]
    │  │     ├─ electrons
    │  │     │  ├─ density ➡ Function [m^-3]
    │  │     │  ├─ density_fast ➡ Function [m^-3]
    │  │     │  ├─ pressure ➡ Function [Pa]
    │  │     │  ├─ pressure_fast_parallel ➡ Function [Pa]
    │  │     │  ├─ pressure_fast_perpendicular ➡ Function [Pa]
    │  │     │  └─ pressure_thermal ➡ Function [Pa]
    │  │     ├─ grid
    │  │     │  ├─ area ➡ Function [m^2]
    │  │     │  ├─ psi ➡ Function [Wb]
    │  │     │  ├─ psi_norm ➡ Function
    │  │     │  ├─ rho_tor_norm ➡ 101-element Vector{Float64}
    │  │     │  ├─ surface ➡ Function [m^2]
    │  │     │  └─ volume ➡ Function [m^3]
    │  │     ├─ j_bootstrap ➡ Function [A/m^2]
    │  │     ├─ j_non_inductive ➡ 101-element Vector{Float64} [A/m^2]
    │  │     ├─ j_ohmic ➡ Function [A/m^2]
    │  │     ├─ j_tor ➡ Function [A/m^2]
    │  │     ├─ j_total ➡ Function [A/m^2]
    │  │     ├─ pressure ➡ Function [Pa]
    │  │     ├─ pressure_ion_total ➡ Function [Pa]
    │  │     ├─ pressure_parallel ➡ Function [Pa]
    │  │     ├─ pressure_perpendicular ➡ Function [Pa]
    │  │     ├─ pressure_thermal ➡ Function [Pa]
    │  │     ├─ t_i_average ➡ Function [eV]
    │  │     └─ time ➡ Function [s]
    │  ├─ time ➡ [2.91] [s]
    │  └─ vacuum_toroidal_field
    │     ├─ b0 ➡ Function [T]
    │     └─ r0 ➡ Function [m]
    ├─ core_sources
    │  └─ vacuum_toroidal_field
    │     ├─ b0 ➡ Function [T]
    │     └─ r0 ➡ Function [m]
    ├─ costing
    │  ├─ cost_decommissioning
    │  │  └─ cost ➡ Function [$M]
    │  ├─ cost_direct_capital
    │  │  └─ cost ➡ Function [$M]
    │  └─ cost_operations
    │     └─ yearly_cost ➡ Function [$M/year]
    ├─ equilibrium
    │  ├─ time ➡ [2.91] [s]
    │  ├─ time_slice
    │  │  └─ 1
    │  │     ├─ boundary
    │  │     │  ├─ elongation ➡ Function
    │  │     │  ├─ elongation_lower ➡ Function
    │  │     │  ├─ elongation_upper ➡ Function
    │  │     │  ├─ geometric_axis
    │  │     │  │  ├─ r ➡ Function [m]
    │  │     │  │  └─ z ➡ Function [m]
    │  │     │  ├─ minor_radius ➡ Function [m]
    │  │     │  ├─ squareness ➡ Function
    │  │     │  ├─ squareness_lower_inner ➡ Function
    │  │     │  ├─ squareness_lower_outer ➡ Function
    │  │     │  ├─ squareness_upper_inner ➡ Function
    │  │     │  ├─ squareness_upper_outer ➡ Function
    │  │     │  ├─ triangularity ➡ Function
    │  │     │  ├─ triangularity_lower ➡ Function
    │  │     │  └─ triangularity_upper ➡ Function
    │  │     ├─ global_quantities
    │  │     │  ├─ energy_mhd ➡ Function [J]
    │  │     │  ├─ ip ➡ 1.24079e+06 [A]
    │  │     │  ├─ magnetic_axis
    │  │     │  │  ├─ b_field_tor ➡ Function [T]
    │  │     │  │  ├─ r ➡ Function [m]
    │  │     │  │  └─ z ➡ Function [m]
    │  │     │  ├─ psi_axis ➡ Function [Wb]
    │  │     │  ├─ psi_boundary ➡ Function [Wb]
    │  │     │  ├─ q_95 ➡ Function
    │  │     │  └─ q_axis ➡ Function
    │  │     ├─ profiles_1d
    │  │     │  ├─ dpressure_dpsi ➡ Function [Pa.Wb^-1]
    │  │     │  ├─ dvolume_drho_tor ➡ 101-element Vector{Float64} [m^2]
    │  │     │  ├─ f ➡ 101-element Vector{Float64} [T.m]
    │  │     │  ├─ geometric_axis
    │  │     │  │  ├─ r ➡ Function [m]
    │  │     │  │  └─ z ➡ Function [m]
    │  │     │  ├─ gm1 ➡ 101-element Vector{Float64} [m^-2]
    │  │     │  ├─ gm2 ➡ 101-element Vector{Float64} [m^-2]
    │  │     │  ├─ gm9 ➡ 101-element Vector{Float64} [m^-1]
    │  │     │  ├─ j_parallel ➡ Function [A/m^2]
    │  │     │  ├─ j_tor ➡ 101-element Vector{Float64} [A.m^-2]
    │  │     │  ├─ psi ➡ 101-element Vector{Float64} [Wb]
    │  │     │  ├─ psi_norm ➡ Function
    │  │     │  ├─ q ➡ 101-element Vector{Float64}
    │  │     │  ├─ rho_tor ➡ 101-element Vector{Float64} [m]
    │  │     │  └─ rho_tor_norm ➡ 101-element Vector{Float64}
    │  │     └─ time ➡ Function [s]
    │  └─ vacuum_toroidal_field
    │     └─ b0 ➡ [-2.08192] [T]
    ├─ stability
    │  └─ all_cleared ➡ Function
    └─ summary
       ├─ fusion
       │  └─ power
       │     └─ value ➡ Function [W]
       ├─ global_quantities
       │  ├─ b0
       │  │  └─ value ➡ Function [T]
       │  ├─ beta_pol_mhd
       │  │  └─ value ➡ Function
       │  ├─ beta_tor
       │  │  └─ value ➡ Function
       │  ├─ beta_tor_mhd
       │  │  └─ value ➡ Function
       │  ├─ beta_tor_norm
       │  │  └─ value ➡ Function
       │  ├─ beta_tor_norm_mhd
       │  │  └─ value ➡ Function
       │  ├─ beta_tor_thermal_norm
       │  │  └─ value ➡ Function
       │  ├─ current_bootstrap
       │  │  └─ value ➡ Function [A]
       │  ├─ current_non_inductive
       │  │  └─ value ➡ Function [A]
       │  ├─ current_ohm
       │  │  └─ value ➡ Function [A]
       │  ├─ energy_thermal
       │  │  └─ value ➡ Function [J]
       │  ├─ h_98
       │  │  └─ value ➡ Function
       │  ├─ ip
       │  │  └─ value ➡ Function [A]
       │  ├─ r0
       │  │  └─ value ➡ Function [m]
       │  ├─ tau_energy
       │  │  └─ value ➡ Function [s]
       │  └─ tau_energy_98
       │     └─ value ➡ Function [s]
       ├─ heating_current_drive
       │  ├─ power_launched_ec
       │  │  └─ value ➡ Function [W]
       │  ├─ power_launched_ic
       │  │  └─ value ➡ Function [W]
       │  ├─ power_launched_lh
       │  │  └─ value ➡ Function [W]
       │  ├─ power_launched_nbi
       │  │  └─ value ➡ Function [W]
       │  └─ power_launched_total
       │     └─ value ➡ Function [W]
       ├─ local
       │  ├─ magnetic_axis
       │  │  ├─ n_e
       │  │  │  └─ value ➡ Function [m^-3]
       │  │  ├─ t_e
       │  │  │  └─ value ➡ Function [eV]
       │  │  ├─ t_i_average
       │  │  │  └─ value ➡ Function [eV]
       │  │  └─ zeff
       │  │     └─ value ➡ Function
       │  └─ separatrix
       │     ├─ n_e
       │     │  └─ value ➡ Function [m^-3]
       │     ├─ t_e
       │     │  └─ value ➡ Function [eV]
       │     ├─ t_i_average
       │     │  └─ value ➡ Function [eV]
       │     └─ zeff
       │        └─ value ➡ Function
       └─ volume_average
          ├─ n_e
          │  └─ value ➡ Function [m^-3]
          ├─ t_e
          │  └─ value ➡ Function [eV]
          ├─ t_i_average
          │  └─ value ➡ Function [eV]
          └─ zeff
             └─ value ➡ Function




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
