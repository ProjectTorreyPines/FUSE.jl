# Transport Solver


```@julia
using Revise
using FUSE
using Plots;
FUSE.logging(Logging.Info);
```

### Let's initialize the case with the DIII-D standard Hmode case

This is a standard DIII-D shot with some 5MW of NBI and some Hmode initialized profiles (not from experiment)


```@julia
dd, ini, act = FUSE.init(:D3D, do_plot=false);
```

### Take a look at the parameters of the actors associated with the transport solver
1. ActorFluxMatcher
2. ActorFluxCalculator
3. ActorTGLF
4. ActorNeoclassical


```@julia
display(act.ActorFluxMatcher)
display(act.ActorFluxCalculator)
display(act.ActorTGLF)
display(act.ActorNeoclassical)
```

#### We can the actors for the transport fluxes individually


```@julia
# We are running with the fast tglfnn model
act.ActorTGLF.sat_rule = :sat0
act.ActorTGLF.electromagnetic = false
act.ActorTGLF.nn = true
FUSE.ActorTGLF(dd,act)
FUSE.ActorNeoclassical(dd,act);
plot(dd.core_transport)
```

#### act.ActorFluxMatcher defines what is evolved. In this case:
   -  Electron Temperature Te
   -  Ion temperature Ti
   -  Electron density ne

We are keeping the rotation fixed and use Deuterium for quasi neutrality and let Carbon match the ne_scale lengths

#### Setting up the actor parameters in act is next


```@julia
# Resetting dd, ini, act so that we don't have to scroll up every time
dd, ini, act = FUSE.init(:D3D, do_plot=false);
#dd, ini, act = FUSE.init(:ITER, init_from=:ods; do_plot=false)

act.ActorTGLF.warn_nn_train_bounds=true
# We are running with the fast tglfnn model
act.ActorTGLF.nn = true
act.ActorTGLF.sat_rule = :sat0
act.ActorTGLF.electromagnetic = false

act.ActorFluxMatcher.rho_transport = 0.3:0.1:0.8
act.ActorFluxMatcher.max_iterations = 100
act.ActorFluxMatcher.optimizer_algorithm = :anderson # or :jacobian_based
act.ActorFluxMatcher.step_size = 0.2
act.ActorFluxMatcher.verbose = true
act.ActorFluxMatcher.evolve_rotation = :fixed

# show pre evolution
display(act.ActorFluxMatcher)
display(plot(dd.core_profiles, label=" before"))

#FUSE.ActorPedestal(dd,act)
actor_transport = FUSE.ActorFluxMatcher(dd, act)

# show after
display(plot!(dd.core_profiles, label=" after"))

# plot the flux_matching 
display(plot(dd.core_transport))
```

### For the channels that we evolved the flux_matching looks spot on!


```@julia
# let's see if our end result satisfies quasi neutral
IMAS.is_quasi_neutral(dd)
```

#### How does TGLF compare to TGLF_nn in this case?


```@julia
act.ActorTGLF.nn = false

# Dialing the iterations a bit down since tglf_sat0 is much slower than it's neural net counterpart
act.ActorFluxMatcher.max_iterations = 30

display(plot(dd.core_profiles, label="  tglfnn"))
FUSE.ActorFluxMatcher(dd, act)

display(plot!(dd.core_profiles, label="  tglf"))
display(plot(dd.core_transport))
```

#### As you can see TGLF sat0 converges well except for a point @ rho = 0.5
#### However, as you can see this doesn't affect the results very much!
