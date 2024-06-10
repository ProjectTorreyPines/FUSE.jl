# FPP multi objective optimization

### Import packages


```@julia
@time using FUSE
using Plots;
gr();
FUSE.logging(Logging.Info; actors=Logging.Error);
```

### Define a working folder


```@julia
save_folder = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain_no0ohm_explore2_minopt_maxflat"
```

### Setup distributed computing environment


```@julia
old_dir = pwd()
mkpath(save_folder)
try
    cd(save_folder) # this is to save temporary distributed files in working folder
    FUSE.parallel_environment("omega", 128 * 2) # ("saga",120)
finally
    cd(old_dir)
end
display(pwd())
using Distributed
@everywhere import FUSE
```

### Get `ini` and `act` for FPP case and custmize as needed


```@julia
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars, STEP=true)
#act.ActorPFcoilsOpt.optimization_scheme = :none; # don't spend time optimizing the PFs
act.ActorStabilityLimits.models = [:beta_troyon_1984, :model_201, :model_401] # include βn check
```

### Define optimization variables and ranges


```@julia
# nominal value and ranges
ini.ec_launchers.power_launched = ini.ec_launchers.power_launched[1] ↔ [1e6, 200e6];
ini.core_profiles.zeff = ini.core_profiles.zeff ↔ [1.1, 2.5]
ini.core_profiles.ne_setting = :greenwald_fraction
ini.core_profiles.ne_value = 1.0 ↔ [0.8, 1.3]
#ini.equilibrium.δ = ini.equilibrium.δ ↔ [-0.7,0.7]
#ini.equilibrium.ζ = ini.equilibrium.ζ ↔ [0,0.2]
ini.equilibrium.κ = 0.95 # k set to be 95% of maximum controllable value
ini.equilibrium.B0 = ini.equilibrium.B0 ↔ [1.0, 20.0]
ini.equilibrium.ip = ini.equilibrium.ip ↔ [1.0e6, 22e6]
ini.equilibrium.R0 = ini.equilibrium.R0 ↔ [ini.equilibrium.R0, 10.0];
ini.requirements.log10_flattop_duration = log10(3600.0) ↔ [log10(3600.0), log10(1000.0 * 3600.0)];
ini.tf.technology = :HTS ↔ (:HTS, :LTS);
ini.oh.technology = :HTS ↔ (:HTS, :LTS);
```

### As a good practice, test the actor/workflow that you want to optimize


```@julia
#dd = FUSE.init(ini, act);
#FUSE.ActorWholeFacility(dd, act);
```


```@julia
#FUSE.digest(dd)
```

### See what are the possible optimization objectives and constraints


```@julia
# FUSE comes with a library of objective and constraints functions
OFL = deepcopy(FUSE.ObjectiveFunctionsLibrary)
CFL = deepcopy(FUSE.ConstraintFunctionsLibrary)
println("== OBJECTIVE FUNCTIONS ==")
display(OFL)
println()
println("== CONSTRAINT FUNCTIONS ==")
display(CFL)
```

## Set the optimization objectives and constraints


```@julia
objective_functions = [OFL[:min_βn], OFL[:min_capital_cost], OFL[:max_log10_flattop]]
constraint_functions = [CFL[:min_required_power_electric_net]]
println("== OBJECTIVE FUNCTIONS ==")
display(objective_functions)
println()
println("== CONSTRAINT FUNCTIONS ==")
display(constraint_functions)
```

### Setup and run optimization


```@julia
# option to resume an optimization where it was left off
if false
    continue_state = state
else
    continue_state = missing
end

# define optimization parameters
# For real optimization studies the population size (N) and number of iterations should be bigger
# eg. N=100, iterations=25
optimization_parameters = Dict(
    :N => max(4, Int(floor((nprocs() - 1) / 2)) * 2), # even number
    :iterations => 100,
    :continue_state => continue_state,
    :save_folder => save_folder)

# run optimization
state = FUSE.workflow_multiobjective_optimization(ini, act, FUSE.ActorWholeFacility, objective_functions, constraint_functions; optimization_parameters...);

```

## Remember to always release computing resources!


```@julia
for i in workers()
    rmprocs(i)
end
```
