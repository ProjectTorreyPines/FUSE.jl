# FPP multi objective optimization

### Setup distributed computing environment

See more details here: https://fuse.help/parallel.html


```julia
if gethostname() == "saga.cluster"
    np = 120
    using Pkg
    Pkg.activate("..")
    using Distributed
    using ClusterManagers
    if nprocs() < np
        addprocs(SlurmManager(np-nprocs()+1),exclusive="", topology=:master_worker)
    end
else
    using Distributed
    np = 4
    if nprocs() < np
        addprocs(np-nprocs()+1, topology=:master_worker)
    end
end
println("Working with $(nprocs()) processes")
```

    Working with 5 processes


### Import packages

NOTE: Import the `plotlyjs()` plotting backend instead of the usual gr() to interactively look at results in 3D


```julia
using Revise
using FUSE
using Plots; gr();
global_logger(FUSE.logger);
```

### Get `ini` and `act` for FPP case and custmize as needed


```julia
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars)
act.ActorTauenn.transport_model = :h98y2
act.ActorPFcoilsOpt.optimization_scheme = :none; # don't spend time optimizing the PFs
```

    WARNING: both ImageMetadata and ImageAxes export "data"; uses of it in module Images must be qualified


### As a good practice, test the actor/workflow that you want to optimize first


```julia
dd = FUSE.init(ini, act)
FUSE.ActorWholeFacility(dd, act);
```


```julia
# define the optimization objectives
OFL = FUSE.ObjectivesFunctionsLibrary
objective_functions = [OFL[:max_fusion], OFL[:min_cost], OFL[:max_flattop]]

# option to resume an optimization where it was left off
if true
    continue_results = missing
else
    continue_results = results
end

# define optimization parameters
# For real optimization studies the population size (N) and number of iterations should be bigger
# eg. N=1000, iterations=100
optimization_parameters = Dict(
    :N => 4,
    :iterations => 10,
    :continue_results => continue_results)

# run optimization
results = FUSE.workflow_multiobjective_optimization(ini, act, FUSE.ActorWholeFacility, objective_functions; optimization_parameters...);
```

    Running on 4 worker processes
    == Actuators ==
    [34m[1mini.ec_launchers.power_launched[22m[39m
    [0m[1m- units: [22m[0mW
    [0m[1m- description: [22m[0mEC launched power
    [0m[1m- value: [22m[0m4.5e7
    [0m[1m- base: [22m[0m4.5e7
    [0m[1m- default: [22m[0mmissing
    [0m[1m- lower: [22m[0m3.0e7
    [0m[1m- upper: [22m[0m1.0e8
    [34m[1mini.core_profiles.zeff[22m[39m
    [0m[1m- units: [22m
    [0m[1m- description: [22m[0mEffective ion charge
    [0m[1m- value: [22m[0m1.1
    [0m[1m- base: [22m[0m1.1
    [0m[1m- default: [22m[0mmissing
    [0m[1m- lower: [22m[0m1.1
    [0m[1m- upper: [22m[0m2.5
    [34m[1mini.core_profiles.greenwald_fraction[22m[39m
    [0m[1m- units: [22m
    [0m[1m- description: [22m[0mGreenwald fraction, ne_vol / ne_gw
    [0m[1m- value: [22m[0m0.9
    [0m[1m- base: [22m[0m0.9
    [0m[1m- default: [22m[0mmissing
    [0m[1m- lower: [22m[0m0.8
    [0m[1m- upper: [22m[0m0.95
    
    == Objectives ==
    [34m[1mmax_fusion[22m[39m [MW]
    [34m[1mmin_cost[22m[39m [$M]
    [34m[1mmax_flattop[22m[39m [hours]
    


    [32mIteration 100%|████████████████████████████| Time: 0:06:23 (38.36  s/it)[39m


    396.822722 seconds (21.55 M allocations: 1.132 GiB, 0.08% gc time, 1.53% compilation time)


### Plot multi-objective optimization results


```julia
plot(results, [1, 2, 3]; design_space=false)
```




    
![svg](multi_objective_optimization_files/multi_objective_optimization_11_0.svg)
    



### How to: Define and use a custom FUSE workflow


```julia
# Here `@everywhere` is needed to make all processes aware of the custom function
@everywhere function workflow_custom(ini,act)
    FUSE.init(dd, ini, act)
    FUSE.ActorEquilibriumTransport(dd, act)
    FUSE.ActorCXbuild(dd, act)
    return dd
end

# results = FUSE.workflow_multiobjective_optimization(ini, act, custom_workflow, objective_functions; optimization_parameters...);
```
