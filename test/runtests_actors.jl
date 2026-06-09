using FUSE
using Test

@testset "fluxmatcher" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    act.ActorFluxMatcher.max_iterations = 2
    act.ActorFluxMatcher.evolve_pedestal = true

    act.ActorFluxMatcher.algorithm = :simple
    act.ActorTGLF.model = :TJLF
    FUSE.ActorFluxMatcher(dd, act)
    act.ActorFluxMatcher.algorithm = :anderson
    act.ActorTGLF.model = :TGLFNN
    FUSE.ActorFluxMatcher(dd, act)
end

@testset "ActorLocking" begin
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    # Small grid and minimal NN for test speed
    N_grid  = 10
    fast_nn = FUSE.NNparams(hidden_sizes=[10], n_epochs=10, batch_size=8)

    # ─── task = :single_case ─────────────────────────────────────────────────
    @testset "single_case" begin
        # Solve one ODE trajectory only — no grid, no NN, no stored results
        actor = FUSE.ActorLocking(dd, act;
                    task          = :single_case,
                    overwrite_params = true)
        @test actor.results === nothing
    end

    # ─── task = :solve_system, control_type = :EF ────────────────────────────
    @testset "solve_system_EF" begin
        actor = FUSE.ActorLocking(dd, act;
                    task             = :solve_system,
                    control_type     = :EF,
                    grid_size        = N_grid,
                    overwrite_params = false,
                    nn_params        = fast_nn)
        r = actor.results

        @test r !== nothing
        # ODE results: one row per grid point, 5 state vars (RP-RW system)
        @test size(r.ode_sols, 1) == N_grid^2
        @test size(r.ode_sols, 2) == 5
        # k-means labels
        @test length(r.locking_labels) == N_grid^2
        @test all(l ∈ (1, 2) for l in r.locking_labels)
        # NN probability model trained and callable
        @test r.prob !== nothing
        C1 = actor.ode_params.Control1[1]
        C2 = actor.ode_params.Control2[1]
        @test 0.0 ≤ r.prob(C1, C2) ≤ 1.0
        # Analytic bifurcation bounds present (NL_saturation_ON = false by default)
        @test r.bifurcation_bounds !== nothing
        @test size(r.bifurcation_bounds) == (N_grid, N_grid)
        @test all(isfinite, r.bifurcation_bounds)
    end

    # ─── task = :solve_system, control_type = :LinStab ───────────────────────
    @testset "solve_system_LinStab" begin
        actor = FUSE.ActorLocking(dd, act;
                    task             = :solve_system,
                    control_type     = :LinStab,
                    grid_size        = N_grid,
                    overwrite_params = true,
                    nn_params        = fast_nn)
        r = actor.results
        @test r !== nothing
        @test size(r.ode_sols, 1) == N_grid^2
        @test r.prob !== nothing
    end

    # ─── task = :solve_system, control_type = :NLsaturation ──────────────────
    @testset "solve_system_NLsaturation" begin
        actor = FUSE.ActorLocking(dd, act;
                    task             = :solve_system,
                    control_type     = :NLsaturation,
                    NL_saturation_ON = true,
                    grid_size        = N_grid,
                    overwrite_params = true,
                    nn_params        = fast_nn)
        r = actor.results
        @test r !== nothing
        # Bifurcation bounds are suppressed when NL saturation is active
        @test r.bifurcation_bounds === nothing
    end

    # ─── task = :calc_prob (retrain NN; ODEs must not be re-solved) ──────────
    @testset "calc_prob_no_ode_rerun" begin
        # Seed the disk file with a known solve
        actor_solve = FUSE.ActorLocking(dd, act;
                          task             = :solve_system,
                          control_type     = :EF,
                          grid_size        = N_grid,
                          overwrite_params = true,
                          nn_params        = fast_nn)
        saved_ode = copy(actor_solve.results.ode_sols)

        # Fresh actor with calc_prob loads the saved ODE results and retrains
        actor_prob = FUSE.ActorLocking(dd, act;
                         task             = :calc_prob,
                         overwrite_params = true,
                         grid_size        = N_grid,
                         nn_params        = fast_nn)
        @test actor_prob.results !== nothing
        @test actor_prob.results.prob !== nothing
        # ODE grid must be identical to what was saved — not re-computed
        @test actor_prob.results.ode_sols ≈ saved_ode
    end

end  # @testset "ActorLocking"


@testset "Mars Actor CHEASE run" begin
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    # control execution of CHEASE & MARS for testing purposes
    act.ActorMars.run_equilibrium = true
    act.ActorMars.run_MHD = false
    act.ActorMars.wall_type = :no_wall # Begin with NO wall
    #act.ActorMars.chease_exec = "/path/to/MarsQ_package/CheaseMerge/chease.x"
    #act.ActorMars.mars_exec = "/path/to/MarsQ_package/MarsQ/marsq.x"

    # Configure CHEASE parameters for testing
    chease_overrides = (NPSI=64, NVEXP=2, NCSCAL=4)

    # configure MARS parameters for testing
    mars_overrides = FUSE.MarsOverrides()
    mars_overrides.BASIC[:M1] = -10
    mars_overrides.BASIC[:NV] = 120 # moves the IW in from the original CHEASE location

    if !isfile(act.ActorMars.chease_exec)
        @test_skip "CHEASE executable not found at $(act.ActorMars.chease_exec). Skipping CHEASE & MARS test."
    else
        # The actor creates and manages its own run directory; keep it across the
        # chained runs below (restart / MHD-only depend on each other's files).
        act.ActorMars.clear_workdir = false

        @info "Test 1: Clean CHEASE equilibrium run"
        actor = FUSE.ActorMars(dd, act; chease_overrides)
        run_dir = actor.par.save_dir
        @testset "CHEASE output files" begin
            @test isfile(joinpath(run_dir, "EXPEQ"))
            @test isfile(joinpath(run_dir, "datain"))
        end

        @info "Test 2: CHEASE with resistive wall & higher beta"
        act.ActorMars.number_surfaces = 2
        act.ActorMars.wall_type = :limiter
        FUSE.ActorMars(dd, act; chease_overrides=(CFBAL=1.1,), save_dir=run_dir)

        @info "Test 3: CHEASE restart from existing equilibrium"
        act.ActorMars.restart_equilibrium = true
        FUSE.ActorMars(dd, act; chease_overrides=(NCSCAL=2,), save_dir=run_dir) # keep Ip fixed

        @info "Test 4: MARS MHD stability run"
        act.ActorMars.run_equilibrium = false # do NOT rerun CHEASE, just run MARS on the existing equilibrium
        act.ActorMars.run_MHD = true
        FUSE.ActorMars(dd, act; mars_overrides, save_dir=run_dir)

        @testset "MARS output files" begin
            @test isfile(joinpath(run_dir, "log_mars"))
            @test isfile(joinpath(run_dir, "RESULT.OUT"))
        end

        @testset "MARS outputs in dd.mhd_linear" begin
            mode = dd.mhd_linear.time_slice[].toroidal_mode[1]

            # growth rate & frequency (converted to SI in _finalize)
            @test mode.n_tor == -1
            @test isfinite(mode.growthrate)
            @test isfinite(mode.frequency)

            # displacement eigenfunction on the (s, m) grid
            ns = length(mode.plasma.grid.dim1)
            nm = length(mode.plasma.grid.dim2)
            @test size(mode.plasma.displacement_perpendicular.real) == (ns, nm)
            @test size(mode.plasma.displacement_perpendicular.imaginary) == (ns, nm)

            # Alfvén-time profile on the radial grid
            @test length(mode.plasma.tau_alfven) == ns
            @test all(mode.plasma.tau_alfven .> 0)

            # real-space R,Z geometry on the (s, χ) grid
            cs = mode.plasma.coordinate_system
            @test size(cs.r) == size(cs.z)
            @test size(cs.r, 1) == length(cs.grid.dim1)
            @test size(cs.r, 2) == length(cs.grid.dim2)

            # plot recipe for the mode structure works
            @test (FUSE.plot(dd.mhd_linear); true)
        end
    end
end
         
