# Poloidal field coil optimization (ITER)


```julia
using Revise
using FUSE
using Plots; gr();
global_logger(FUSE.logger);
```

### Initialization of `dd`, `ini`, `act` for the ITER use case


```julia
dd, ini, act = FUSE.init(:ITER, init_from=:ods, do_plot=true);
```

    WARNING: both ImageMetadata and ImageAxes export "data"; uses of it in module Images must be qualified
    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mdd.dataset_description was skipped in IMAS data dictionary
    [33m[1m└ [22m[39m[90m@ IMASDD ~/.julia/dev/IMASDD/src/data.jl:1067[39m



    
![svg](PF_coil_optimization_files/PF_coil_optimization_3_1.svg)
    



    
![svg](PF_coil_optimization_files/PF_coil_optimization_3_2.svg)
    



    [1m13×8 DataFrame[0m
    [1m Row [0m│[1m group  [0m[1m name          [0m[1m ΔR      [0m[1m R_start [0m[1m R_end   [0m[1m material  [0m[1m area     [0m[1m volume    [0m
    [1m     [0m│[90m String [0m[90m String        [0m[90m Float64 [0m[90m Float64 [0m[90m Float64 [0m[90m String    [0m[90m Float64  [0m[90m Float64   [0m
    ─────┼──────────────────────────────────────────────────────────────────────────────────
       1 │ in                        0.8      0.0      0.8              10.7597     27.0421
       2 │ in      OH                1.3      0.8      2.1   Nb3Sn      17.4845    159.295
       3 │ hfs     TF                1.1      2.1      3.2   Nb3Sn      36.7172    484.472
       4 │ hfs     vacuum vessel     0.3      3.2      3.5   Water      20.8103    934.786
       5 │ hfs     shield            0.4      3.5      3.9   Tungsten    9.81998   356.626
       6 │ hfs     wall              0.06     3.9      3.96  Steel       4.17208   160.693
       7 │ lhfs    plasma            4.4      3.96     8.36  DT_plasma  28.8365   1062.79
       8 │ lfs     wall              0.17     8.36     8.53  Steel       4.17208   160.693
       9 │ lfs     shield            0.4      8.53     8.93  Tungsten    9.81998   356.626
      10 │ lfs     vacuum vessel     1.05     8.93     9.98  Water      20.8103    934.786
      11 │ lfs     TF                1.1      9.98    11.08  Nb3Sn      36.7172    484.472
      12 │ out                       2.34    11.08    13.42             92.1117   4795.84
      13 │ out     cryostat          0.3     13.42    13.72  Steel       8.89184   566.361



    
![svg](PF_coil_optimization_files/PF_coil_optimization_3_4.svg)
    



    
![svg](PF_coil_optimization_files/PF_coil_optimization_3_5.svg)
    



    
![svg](PF_coil_optimization_files/PF_coil_optimization_3_6.svg)
    



    
![svg](PF_coil_optimization_files/PF_coil_optimization_3_7.svg)
    


### Run the SteadyStateCurrent actor in order to estimate how much ohmic current will be required during flattop
[SteadyStateCurrent actor](https://fuse.help/actors.html#SteadyStateCurrent)



```julia
FUSE.ActorSteadyStateCurrent(dd,act);
```

### Estimate how much flux is required during start-up
[FluxSwing actor](https://fuse.help/actors.html#FluxSwing)




```julia
FUSE.ActorFluxSwing(dd,act; operate_at_j_crit=true,j_tolerance=0.0)
# critial_j shoudl be the same as max_j
dd.build.oh
# Estimated flattop durration at maximum oh current operation suggests that ITER can run for 3500 seconds! 
```




    [0m[1moh[22m
    ├─ [0mcritical_b_field[31m ➡ [39m[31m18.7235[39m
    ├─ [0mcritical_j[31m ➡ [39m[31m9.24932e+06[39m
    ├─ [0mflattop_duration[31m ➡ [39m[33m1000[39m
    ├─ [0mflattop_estimate[31m ➡ [39m[31m3466.12[39m
    ├─ [0mmax_b_field[31m ➡ [39m[31m15.0455[39m
    ├─ [0mmax_j[31m ➡ [39m[31m9.2099e+06[39m
    └─ [0m[1mtechnology[22m
       ├─ [0mJxB_strain[31m ➡ [39m[31m-0.05[39m
       ├─ [0mfraction_stainless[31m ➡ [39m[31m0.46[39m
       ├─ [0mfraction_void[31m ➡ [39m[31m0.1[39m
       ├─ [0mmaterial[31m ➡ [39m[35m"Nb3Sn"[39m
       ├─ [0mratio_SC_to_copper[31m ➡ [39m[31m1[39m
       ├─ [0mtemperature[31m ➡ [39m[31m4.2[39m
       └─ [0mthermal_strain[31m ➡ [39m[31m-0.64[39m




### Run PF coils optimization and plot


```julia
FUSE.init_pf_active(dd,ini,act)
#plot(dd.pf_active)
#act.ActorPFcoilsOpt[:optimization_scheme]
FUSE.ActorPFcoilsOpt(dd,act;optimization_scheme=:currents,do_plot=true);
```


    
![svg](PF_coil_optimization_files/PF_coil_optimization_9_0.svg)
    



    
![svg](PF_coil_optimization_files/PF_coil_optimization_9_1.svg)
    



    
![svg](PF_coil_optimization_files/PF_coil_optimization_9_2.svg)
    



    
![svg](PF_coil_optimization_files/PF_coil_optimization_9_3.svg)
    

