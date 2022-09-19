# Validation of the coupled transport-equilibrium workflow


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Running one validation sample per tokamak in the HDB5 dataset


```@julia
FUSE.workflow_HDB5_validation(n_samples_per_tokamak=1, show_dd_plots=true, plot_database=true, show_progress=false);
```
