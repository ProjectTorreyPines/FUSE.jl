# FPP multi objective optimization

### Import packages


```@julia
using Revise
using FUSE
using Plots; gr();
FUSE.logging(Logging.Info; actors=Logging.Error);
```

### Setup distributed computing environment

See more details here: https://fuse.help/parallel.html


```@julia
FUSE.parallel_environment("localhost")
using Distributed
```

### Get `ini` and `act` for FPP case and custmize as needed


```@julia
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars)
act.ActorPFcoilsOpt.optimization_scheme = :none; # don't spend time optimizing the PFs
```

### As a good practice, test the actor/workflow that you want to optimize first


```@julia
dd = FUSE.init(ini, act)
FUSE.ActorWholeFacility(dd, act);
```


```@julia
#IMAS.freeze(dd.balance_of_plant)
plot(dd.core_profiles)
```

### Define optimization variables and ranges


```@julia
# nominal value and ranges
ini_opt = deepcopy(ini)
ini_opt.ec_launchers.power_launched = ini.ec_launchers.power_launched ↔ [10e6, 200e6];
ini_opt.core_profiles.zeff = ini.core_profiles.zeff ↔ [1.1, 2.5]
ini_opt.equilibrium.δ = ini.equilibrium.δ ↔ [-0.7,0.7]
ini_opt.equilibrium.ζ = ini.equilibrium.ζ ↔ [0,0.2]
ini_opt.equilibrium.κ = ini.equilibrium.κ ↔ [1.3,2.5]
ini_opt.equilibrium.B0 = ini.equilibrium.B0 ↔ [1.0, 20.]
ini_opt.equilibrium.ip = ini.equilibrium.ip ↔ [1.0e6, 22e6]
ini_opt.equilibrium.R0 = ini.equilibrium.R0 ↔ [ini.equilibrium.R0, 10.0];
```

### Define the optimization objectives


```@julia
# FUSE comes with a library of objective functions
OFL = deepcopy(FUSE.ObjectivesFunctionsLibrary)
OFL[:max_power_electric_net].target = 200.0
objective_functions = [OFL[:max_power_electric_net], OFL[:min_log10_levelized_CoE], OFL[:max_log10_flattop]]

CFL = deepcopy(FUSE.ConstraintFunctionsLibrary)
CFL[:target_power_electric_net].limit = 200.0
CFL[:target_power_electric_net].tolerance = 0.01
constraint_functions = [CFL[:target_power_electric_net], CFL[:steady_state]]

# ...but one can define custom objectives and constraints too
# target_power_electric = FUSE.ObjectiveFunction(:target_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net)/1E6, 200)
#objective_functions = [target_power_electric, OFL[:min_cost], OFL[:max_flattop]]

display(objective_functions)
display(constraint_functions)
```

### Setup and run optimization


```@julia
# option to resume an optimization where it was left off
if false
    continue_results = results
else
    continue_results = missing
end

# define optimization parameters
# For real optimization studies the population size (N) and number of iterations should be bigger
# eg. N=100, iterations=25
optimization_parameters = Dict(
    :N => max(4, Int(floor((nprocs()-1)/2))*2), # even number
    :iterations => 10,
    :continue_results => continue_results,
    :save_folder => "optimization_runs")

# run optimization
results = FUSE.workflow_multiobjective_optimization(ini_opt, act, FUSE.ActorWholeFacility, objective_functions, constraint_functions; optimization_parameters...);
```

### Save optimization results to file


```@julia
# Optimization results can be re-loaded this way:
filename = "optimization_runs/optimization.bson"
@time results = FUSE.load_optimization(filename);
```

### Plot multi-objective optimization results


```@julia
extract_inputs = Dict(
    :Paux => (dd, ini, act) -> ini.ec_launchers.power_launched / 1E6,
    :zeff => (dd, ini, act) -> ini.core_profiles.zeff,
    :κ => (dd, ini, act) -> ini.equilibrium.κ,
    :δ => (dd, ini, act) -> ini.equilibrium.δ,
    :ζ => (dd, ini, act) -> ini.equilibrium.ζ,
    :B0 => (dd, ini, act) -> ini.equilibrium.B0,
    :ip => (dd, ini, act) -> ini.equilibrium.ip / 1E6,
    :R0 => (dd, ini, act) -> ini.equilibrium.R0,
)

extract_outputs = Dict(
    :beta_n => (dd, ini, act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
    :log10_levelized_CoE => (dd, ini, act) -> log10(dd.costing.levelized_CoE),
    :Pfusion => (dd, ini, act) -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6,
    :Pelectric => (dd, ini, act) -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6,
    :log10_flattop => (dd, ini, act) -> log10(dd.build.oh.flattop_duration / 3600.0)
)

path="optimization_runs"
dirs = filter(isdir,readdir(path; join=true));
inputs,outputs = FUSE.load(dirs, [extract_inputs, extract_outputs]);
```


```@julia
using Plots
display(plot([histogram(inputs[:,name], title=name, label="") for name in names(inputs)]...))
plot([histogram(outputs[:,name], title=name, label="") for name in names(outputs)]...)
```


```@julia
xname="log10_levelized_CoE"
x=outputs[:,xname]
yname="Pelectric"
y=outputs[:,yname]
cname="beta_n"
c=outputs[:,cname]
scatter(x, y, marker_z=c ,xlabel=xname, ylabel=yname, colorbar_title=cname, marker=:circle, markerstrokewidth=0, label="")
```

### How to: Define and use a custom FUSE workflow


```@julia
# Here `@everywhere` is needed to make all processes aware of the custom function
@everywhere function workflow_custom(ini, act)
    FUSE.init(dd, ini, act)
    FUSE.ActorEquilibriumTransport(dd, act)
    FUSE.ActorCXbuild(dd, act)
    return dd
end

# results = FUSE.workflow_multiobjective_optimization(ini, act, custom_workflow, objective_functions; optimization_parameters...);
```
