# Neutronics example


```julia
using Revise
using FUSE
using Plots; gr();
global_logger(FUSE.logger);
```

### Initialize FPP v1_demount case
[FPP v1 demount case documentation](https://fuse.help/cases.html#FPP)


```julia
dd, ini, act = FUSE.init(:FPP, version=:v1_demount, init_from=:scalars, do_plot=false);
```

    WARNING: both ImageMetadata and ImageAxes export "data"; uses of it in module Images must be qualified


### Run the Neutronics actor


```julia
FUSE.ActorNeutronics(dd, act; do_plot=true);
```


    
![svg](neutronics_files/neutronics_5_0.svg)
    


### Plot neutron wall loading


```julia
plot(dd.neutronics.time_slice[].wall_loading; cx=true)
```




    
![svg](neutronics_files/neutronics_7_0.svg)
    




```julia
plot(dd.neutronics.time_slice[].wall_loading; cx=false)
```




    
![svg](neutronics_files/neutronics_8_0.svg)
    


