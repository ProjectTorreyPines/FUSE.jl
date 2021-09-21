# FUSE

Getting started
===============

1. Bring in Project Torrey Pines private repos in your Julia development directory:
    ```julia
    using Pkg
    Pkg.develop(url="git@github.com:ProjectTorreyPines/FUSE.jl.git")
    Pkg.develop(url="git@github.com:ProjectTorreyPines/IMAS.jl.git")
    Pkg.develop(url="git@github.com:ProjectTorreyPines/Equilibrium.jl.git")
    ```
    this will clone the repos in the the `FUSE`, `IMAS`, `Equilibrium` development folders under `~/.julia/dev`
  
2. Activate the FUSE environment:
    ```julia
    using Pkg
    Pkg.activate("..")
    using Revise # useful for development
    using FUSE
    using FUSE.IMAS
    ```
3. Play:

    ```julia
    R0 = 1.8
    δ = 0.5
    ϵ = 0.3
    κ = 1.9
    B0 = 2.0
    qstar = 1.5
    beta_t = 0.01

    # initialize equilibrium IDS
    eq0 = FUSE.init(IMAS.equilibrium(), 0.0; B0, R0, ϵ, δ, κ, beta_t, qstar)
    
    # instantiate equilibrium actor
    eqactor = FUSE.SolovevEquilibriumActor(eq0, 0.0)
    
    # step (aka run, solve) your actor
    FUSE.step(eqactor; abs_error=1E-4, verbose=false)
    
    # translate actor internals to equilibrium IDS
    eq1 = FUSE.finalize(eqactor,129)
    ```
