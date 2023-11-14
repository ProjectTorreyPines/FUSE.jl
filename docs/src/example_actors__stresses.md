# Stresses


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Initialize FPP v1_demount case
[FPP v1 demount case documentation](https://fuse.help/cases.html#FPP)


```@julia
ini, act = FUSE.case_parameters(:FPP, version=:v1_demount, init_from=:scalars);
```


```@julia
# FPP TF and OH are free-standing
ini.center_stack
```


```@julia
# Here is the radial build
ini.build.layers
```


```@julia
# remove gap between TF and OH, so that we can buck them together
pop!(ini.build.layers, "gap_TF_OH");
ini.build.layers
```


```@julia
# initialize
dd = FUSE.init(ini, act);
```


```@julia
# print the radial build as a table
display(dd.build.layer)
# plot the radial build (OH and hfs TF are touching)
plot(dd.build; cx=false)
```

### Run some of the actors needed to evaluate the stresses


```@julia
FUSE.ActorEquilibriumTransport(dd, act)
FUSE.ActorFluxSwing(dd, act);
```

### Run the Stresses actor
[ActorStresses documentation](https://fuse.help/actors.html#Stresses)


```@julia
act.ActorStresses.n_points = 51

dd.solid_mechanics.center_stack.plug = 0
dd.solid_mechanics.center_stack.bucked = 0
dd.solid_mechanics.center_stack.noslip = 0

plot(size=(800, 500))

# free standing
FUSE.ActorStresses(dd, act)
plot!(dd.solid_mechanics.center_stack.stress, color=:blue)

# add bucking
dd.solid_mechanics.center_stack.bucked = 1
FUSE.ActorStresses(dd, act)
plot!(dd.solid_mechanics.center_stack.stress, color=:red)

# add plug
dd.solid_mechanics.center_stack.plug = 1
FUSE.ActorStresses(dd, act)
plot!(dd.solid_mechanics.center_stack.stress, color=:green)
```
