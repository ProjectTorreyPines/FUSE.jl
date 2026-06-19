using FUSE
using Test

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
                    overwrite_params = false)
        r = actor.results

        @test r !== nothing
        # ODE results: one row per grid point, 5 state vars (RP-RW system)
        @test size(r.ode_sols, 1) == N_grid^2
        @test size(r.ode_sols, 2) == 5
        # k-means labels
        @test length(r.locking_labels) == N_grid^2
        @test all(l ∈ (1, 2) for l in r.locking_labels)
        # NN training is deferred to task=:calc_prob — :solve_system leaves prob unset
        @test r.prob === nothing
        # Analytic bifurcation bounds present (NL_saturation = false by default)
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
                    ode_params       = (Control2_min=-3.0, Control2_max=-1.0))
        r = actor.results
        @test r !== nothing
        @test size(r.ode_sols, 1) == N_grid^2
        # NN training is deferred to task=:calc_prob — :solve_system leaves prob unset
        @test r.prob === nothing
    end

    # ─── task = :solve_system, control_type = :NLsaturation ──────────────────
    @testset "solve_system_NLsaturation" begin
        actor = FUSE.ActorLocking(dd, act;
                    task             = :solve_system,
                    control_type     = :NLsaturation,
                    NL_saturation    = true,
                    grid_size        = N_grid,
                    overwrite_params = true,
                    ode_params       = (Control2_min=0.05, Control2_max=0.5))
        r = actor.results
        @test r !== nothing
        
    end

    # ─── task = :calc_prob (retrain NN; ODEs must not be re-solved) ──────────
    @testset "calc_prob_no_ode_rerun" begin
        # Seed the disk file with a known solve
        actor_solve = FUSE.ActorLocking(dd, act;
                          task             = :solve_system,
                          control_type     = :EF,
                          grid_size        = N_grid,
                          overwrite_params = true)
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

    # ─── task = :transfer_learning ───────────────────────────────────────────
    @testset "transfer_learning" begin
        # Seed a base NN model on disk via a full :solve_system run
        actor_base = FUSE.ActorLocking(dd, act;
                         task             = :solve_system,
                         control_type     = :EF,
                         grid_size        = N_grid,
                         overwrite_params = true)
        # NN training is deferred to task=:calc_prob — train explicitly before saving
        FUSE.train_locking_nn(actor_base, fast_nn)
        FUSE.save_locking_nn(actor_base)

        # Fine-tune (last layer only) on a focused/sparser control-space sweep
        tl_nn = FUSE.NNparams(hidden_sizes=[10], n_epochs=5, batch_size=8)
        actor_tl = FUSE.ActorLocking(dd, act;
                       task             = :transfer_learning,
                       control_type     = :EF,
                       grid_size        = 5,
                       overwrite_params = true,
                       nn_params        = tl_nn,
                       ode_params       = (Control1_min=1.0, Control1_max=3.0))
        r = actor_tl.results

        @test r !== nothing
        @test size(r.ode_sols, 1) == 5^2
        @test r.prob isa FUSE.LockingNNModel
        C1 = actor_tl.ode_params.Control1[1]
        C2 = actor_tl.ode_params.Control2[1]
        @test 0.0 ≤ r.prob(C1, C2) ≤ 1.0
        # Focused Control1 (Ω0) range was honored
        @test extrema(actor_tl.ode_params.Control1) == (1.0, 3.0)
    end

end  # @testset "ActorLocking"
