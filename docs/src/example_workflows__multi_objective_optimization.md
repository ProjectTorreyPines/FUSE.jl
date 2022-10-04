# FPP multi objective optimization

### Setup distributed computing environment

See more details here: https://fuse.help/parallel.html


```@julia
if gethostname() == "saga.cluster"
    nodes = 4
    np = 30 * nodes
    using Pkg
    Pkg.activate("..")
    using Distributed
    using ClusterManagers
    ENV["JULIA_WORKER_TIMEOUT"] = "180"
    if nprocs() < np
        addprocs(SlurmManager(np - nprocs()), exclusive="", topology=:master_worker)
    end
else
    using Distributed
    np = 4
    if nprocs() < np + 1
        addprocs(np - nprocs() + 1, topology=:master_worker)
    end
end
println("Working with $(nprocs()) processes")
```

### Import packages


```@julia
using Revise
using FUSE
using Plots;
gr();
FUSE.logging(Logging.Info; actors=Logging.Error);
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
ini.core_profiles.zeff = 2.0 ↔ [1.1, 2.5]
ini.core_profiles.greenwald_fraction = 0.9 ↔ [0.8, 0.95]
ini.ec_launchers.power_launched = 45e6 ↔ [20e6, 60e6];
```

### Define the optimization objectives


```@julia
# FUSE comes with a library of objective functions
OFL = FUSE.ObjectivesFunctionsLibrary
objective_functions = [OFL[:max_power_electric_net], OFL[:min_cost], OFL[:max_log10_flattop]]

# ...but one can define custom optimization objectives too
# target_power_electric = FUSE.ObjectiveFunction(:target_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net)/1E6, 200)
# objective_functions = [target_power_electric, OFL[:min_cost], OFL[:max_flattop]]
```

### Setup and run optimization


```@julia
# option to resume an optimization where it was left off
if true
    continue_results = missing
else
    continue_results = results
end

# define optimization parameters
# For real optimization studies the population size (N) and number of iterations should be bigger
# eg. N=100, iterations=25
optimization_parameters = Dict(
    :N => 4, # even number
    :iterations => 10,
    :continue_results => continue_results)

# run optimization
#results = FUSE.workflow_multiobjective_optimization(ini, act, FUSE.ActorWholeFacility, objective_functions; optimization_parameters...);
results = FUSE.workflow_multiobjective_optimization(ini, act, Val{:remote}, objective_functions; optimization_parameters...);
```

### Save optimization results to file


```@julia
filename = "optimization.bson"
#@time FUSE.save_optimization(filename, results)

# Optimization results can be re-loaded this way:
#@time results = FUSE.load_optimization(filename);
```

### Plot multi-objective optimization results


```@julia
using Plots;
gr();
using Interact

design_space = false
pareto = false
if design_space
    xlim = [results.opt_ini[1].lower, results.opt_ini[1].upper]
    ylim = [results.opt_ini[2].lower, results.opt_ini[2].upper]
    zlim = [results.opt_ini[3].lower, results.opt_ini[3].upper]
else
    xlim = [0, 200]
    ylim = [1800, 1900]
    zlim = [0, 5]
end
@manipulate for iteration in 1:25
    iterations = iteration:iteration
    p = []
    for k in [1, 2, 3]
        push!(p, plot(results, [1, 2, 3]; design_space, pareto, color_by=k, max_samples=nothing, iterations, xlim, ylim, zlim, labelfontsize=8, titlefontsize=10, margin=5Plots.mm))
    end
    plot(p..., layout=(1, 3), size=(1200, 300))
end

```


```@julia
# Import the `plotlyjs()` plotting backend instead of the usual `gr()`to interactively look at results in 3D
#using Plots; plotlyjs();
#pareto=true
#plot(results, [1,2,3]; design_space, pareto, color_by=2, max_samples=nothing, iterations=25:25)
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
