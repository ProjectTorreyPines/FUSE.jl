# FUSE

Getting started
===============

1. Bring in Project Torrey Pines private repos in your Julia development directory:
    ```julia
    ] develop "git@github.com:ProjectTorreyPines/FUSE.jl.git"
    ] develop "git@github.com:ProjectTorreyPines/IMAS.jl.git"
    ] develop "git@github.com:ProjectTorreyPines/Equilibrium.jl.git"
    ] develop "git@github.com:ProjectTorreyPines/AD_GS.jl.git"
    ```
    this will clone the repos in the the `FUSE`, `IMAS`, `Equilibrium` development folders under `~/.julia/dev`
  
2. Activate the FUSE environment (NOTE: if not starting from scratch, it may be necessary to remove the `Manifest.toml` file from the development repositories):
    ```bash
    cd ~/.julia/dev/FUSE
    ```
    
    ```julia
    ] activate .    # activate the ~/.julia/dev/FUSE/Project.toml environment
    ] status        # see where we are
    # We need to tell the environment not to look for our packages in our dev folder and not in the official Julia package registry
    ] develop Equilibrium IMAS AD_GS
    ] resolve       # Download packages, update
    ```

4. Play:

    ```julia
    using Pkg
    Pkg.activate(ENV["HOME"] * "/.julia/dev/FUSE")
    using Revise
    using FUSE
    using FUSE.IMAS
    using Plots
    
    R0 = 1.8
    δ = 0.5
    ϵ = 0.3
    κ = 1.9
    B0 = 2.0
    qstar = 1.5
    beta_n = 1.5
    ip = 1e6
    x_point = false

    # initialize one time-slice of the equilibrium IDS
    dd = IMAS.dd()
    resize!(dd.equilibrium.time_slice,1)
    FUSE.init(dd.equilibrium.time_slice[1]; B0, R0, ϵ, δ, κ, beta_n, ip, x_point)

    # instantiate equilibrium actor
    eqactor = FUSE.SolovevEquilibriumActor(eq0)
    
    # step (a.k.a. run, solve) your actor
    FUSE.step(eqactor)
    
    # translate actor internals to equilibrium IDS
    eq1 = FUSE.finalize(eqactor)
    
    # plot equilibrium
    plot(eq1)
    ```
