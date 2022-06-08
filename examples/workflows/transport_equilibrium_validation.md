# Validation of the coupled transport-equilibrium workflow


```julia
using Revise
using FUSE
using Plots; gr();
global_logger(FUSE.logger);
```

### Running 2 samples of the validation workflow of the HDB5 dataset


```julia
FUSE.workflow_HDB5_validation(n_samples_per_tokamak=1,show_dd_plots=true, plot_database=true,show_progress=false);
```

    WARNING: both ImageMetadata and ImageAxes export "data"; uses of it in module Images must be qualified



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_1.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_2.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_3.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_4.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_5.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_6.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_7.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_8.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_9.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_10.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_11.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_12.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_13.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_14.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_15.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_16.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_17.svg)
    



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_18.svg)
    


    Failed runs: 1 out of 6
    Mean Relative error 38.09%



    
![svg](transport_equilibrium_validation_files/transport_equilibrium_validation_3_20.svg)
    


    R² = 0.69, mean_relative_error = 38.09)


### Comment on the results
