# FPP parametric scan

## Import packages


```@julia
using Revise
using FUSE
using Plots; gr();
FUSE.logging(Logging.Info; actors=Logging.Error);
```

    [ Info: Precompiling FUSE [e64856f0-3bb8-4376-b4b7-c03396503992]


### Setup distributed computing environment

See more details here: https://fuse.help/parallel.html


```@julia
FUSE.parallel_environment("localhost")
using Distributed
@everywhere using FUSE
```

### Define workflow and extent of the scan


```@julia
@everywhere function workflow(B0, spf)::IMAS.dd
    dd = IMAS.dd()
    ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars);
    ini.equilibrium.ip *= B0/ini.equilibrium.B0
    ini.equilibrium.B0 = B0
    ini.core_profiles.polarized_fuel_fraction = spf
    try
        FUSE.init(dd, ini, act);
        FUSE.ActorWholeFacility(dd, act);
        return dd
    catch e
        @warn e
        return dd
    end
end

ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars);
B0_ = LinRange(5,6.5,8)
spf_ = LinRange(0,1,10)
cases = [[B0,spf] for B0 in B0_, spf in spf_]
```

### Run simulations


```@julia
using ProgressMeter
DD = @showprogress pmap(case -> workflow(case...), cases);
```

### Extract quantities of interest from series of `dd`'s


```@julia
results=FUSE.extract(reshape(DD,length(DD)); filter_invalid=false)
```

### Plot scan results


```@julia
results_ = reshape(results.Pfusion,(length(spf_),length(B0_)))
scatter(spf_,B0_,results_,xlabel="SPF fraction",ylabel="B0",levels=10,grid=false)
vline!(spf_,color=:gray,alpha=0.2,label="")
hline!(B0_,color=:gray,alpha=0.2,label="")
```




    
![svg](assets/parametric_scan_files/parametric_scan_12_0.svg)
    




```@julia
p=plot(DD[1].equilibrium)
for k in 6:10
    plot!(DD[k].equilibrium)
end
p
```




    
![svg](assets/parametric_scan_files/parametric_scan_13_0.svg)
    




```@julia

```
