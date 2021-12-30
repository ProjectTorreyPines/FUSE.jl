# FUSE.jl

Getting started
===============

1. Clone the FUSE repository in your Julia development directory:
    ```julia
    ] develop "git@github.com:ProjectTorreyPines/FUSE.jl.git"
    ```
    this will clone FUSE.jl in the the `FUSE` development folders under `~/.julia/dev`
  
2. Install all FUSE dependencies:
    ```bash
    cd ~/.julia/dev/FUSE
    make install
    ```

3. Update FUSE and dependencies to latest version
    ```bash
    cd ~/.julia/dev/FUSE
    make update
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
