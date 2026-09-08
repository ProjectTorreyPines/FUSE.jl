using FUSE
using Test

@testset verbose = true "ActorLocking" begin
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    # UPDATE sources to chosen time,if fetching a shot from MDS
    #dd.global_time=3.0
    #act.ActorSources.nb_model=:NBsimple
    #FUSE.ActorSources(dd, act)  

    # Small grid and minimal NN for test speed
    N_grid  = 10
    fast_nn = FUSE.NNparams(hidden_sizes=[10], n_epochs=10, batch_size=8)

    # ─── task = :single_case ─────────────────────────────────────────────────
    @testset "single_case, one ODE trajectory" begin
        # Solve one ODE trajectory only — no grid, no NN, no stored results
        actor = FUSE.ActorLocking(dd, act;
                    task          = :single_case,
                    overwrite_params = true)
        @test actor.results === nothing
    end

    # ─── task = :solve_system, control_type = :EF ────────────────────────────
    @testset "solve_system, error field (EF) control" begin
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
    @testset "solve_system, stability index (LinStab) control" begin
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
    @testset "solve_system, saturation (NLsaturation) control" begin
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
    @testset "calc_prob, retrain NN without re-solving" begin
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

    # ─── task = :eval_prob ───────────────────────────────────────────────────
    @testset "eval_prob, probability at chosen times" begin
        # Seed disk with ODE results (:solve_system) and a probability model
        # (:calc_prob) — :eval_prob loads both rather than recomputing
        FUSE.ActorLocking(dd, act;
            task             = :solve_system,
            control_type     = :EF,
            grid_size        = N_grid,
            overwrite_params = true)
        FUSE.ActorLocking(dd, act;
            task             = :calc_prob,
            control_type     = :EF,
            grid_size        = N_grid,
            overwrite_params = true,
            nn_params        = fast_nn)

        # op_times must be ≥ dd.global_time: earlier testsets already appended an
        # mhd_linear slice there, and resize! rejects a time below the last one
        t_ref    = dd.global_time
        op_times = t_ref .+ [0.0, 1e-3, 2e-3]

        # A scalar op_C2 broadcasts across every op_time
        actor = FUSE.ActorLocking(dd, act;
                    task             = :eval_prob,
                    control_type     = :EF,
                    grid_size        = N_grid,
                    overwrite_params = true,
                    op_times         = op_times,
                    op_C2            = 5.0)

        ev = actor.eval
        @test ev !== nothing
        @test length(ev.times) == length(op_times)
        @test length(ev.C1)    == length(op_times)
        @test ev.C2 == fill(5.0, length(op_times))   # scalar broadcast
        @test all(isfinite, ev.C1)

        # One mhd_linear time_slice per op time, each carrying P(locked) —
        # previously only the last time survived
        for t in op_times
            dd.global_time = t
            mode = dd.mhd_linear.time_slice[].toroidal_mode[1]
            @test 0.0 ≤ mode.stability_metric ≤ 1.0
        end
        dd.global_time = t_ref

        # Probability limit written over the same time base
        probs = filter(m -> occursin("probability", m.identifier.name), dd.limits.model)
        @test length(probs) == 1
        @test length(probs[1].fraction) == length(op_times)
        @test all(isfinite, probs[1].fraction)

        # A vector op_C2 must match op_times one-for-one
        @test_throws ErrorException FUSE.ActorLocking(dd, act;
            task         = :eval_prob,
            control_type = :EF,
            grid_size    = N_grid,
            op_times     = op_times,
            op_C2        = [1.0, 5.0])
    end

    # ─── task = :transfer_learning ───────────────────────────────────────────
    @testset "transfer_learning, fine-tune on focused grid" begin
        # Seed a base NN model on disk via a full :solve_system run
        actor_base = FUSE.ActorLocking(dd, act;
                         task             = :solve_system,
                         control_type     = :EF,
                         grid_size        = N_grid,
                         overwrite_params = true,
                         nn_params        = fast_nn)
        # NN training is deferred to task=:calc_prob — train explicitly before saving
        FUSE.train_locking_nn(actor_base)
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


    # ─── control_type = :LinStab — br (Gauss) → Δt inversion ─────────────────
    @testset "eval_prob, br to Δt inversion (LinStab)" begin
        # Control2_min/max are Δ_RW for :LinStab and must be < 0; the ODEparams
        # defaults are positive, so they have to be given explicitly
        drw_range = (Control2_min=-3.5, Control2_max=-0.05)

        FUSE.ActorLocking(dd, act;
            task             = :solve_system,
            control_type     = :LinStab,
            grid_size        = N_grid,
            overwrite_params = true,
            ode_params       = drw_range)
        FUSE.ActorLocking(dd, act;
            task             = :calc_prob,
            control_type     = :LinStab,
            grid_size        = N_grid,
            overwrite_params = true,
            ode_params       = drw_range,
            nn_params        = fast_nn)

        # later than the times the :EF eval_prob testset already appended
        t_ref    = dd.global_time
        op_times = t_ref .+ [3e-3, 4e-3]

        # op_C2 is the n=1 br amplitude in Gauss here, one per time
        actor = FUSE.ActorLocking(dd, act;
                    task             = :eval_prob,
                    control_type     = :LinStab,
                    grid_size        = N_grid,
                    overwrite_params = true,
                    ode_params       = drw_range,
                    op_times         = op_times,
                    op_C2            = [10.0, 20.0])

        ev = actor.eval
        @test ev !== nothing
        @test ev.C2 == [10.0, 20.0]            # stored in user units, not Δt

        op = actor.ode_params
        # error_field is derived from par.error_field (Gauss), never accumulated
        @test op.error_field ≈ act.ActorLocking.error_field * 1e-4 / actor.par.b0 *
                               op.control_surf / actor.par.m_pol

        # br → Δt must stay below the marginal Δt_crit, and rise with br:
        # larger br ⇒ smaller |Δ_RW| ⇒ Δt closer to Δt_crit
        Δt_crit = op.l21 * op.l12 / op.DeltaW
        Δt = [FUSE._op_C2_to_control2(actor, c2) for c2 in ev.C2]
        @test all(Δt .< Δt_crit)
        @test Δt[2] > Δt[1]

        # Δ_RW = Δt - Δt_crit must be strictly negative (weak stability)
        @test all(Δt .- Δt_crit .< 0.0)

        # one mhd_linear slice per op time, as for :EF
        for t in op_times
            dd.global_time = t
            @test 0.0 ≤ dd.mhd_linear.time_slice[].toroidal_mode[1].stability_metric ≤ 1.0
        end
        dd.global_time = t_ref
    end

    # ─── prob_method engine selection and window validation ──────────────────
    @testset "calc_prob, prob_method engines" begin
        # an even window has no centre point — rejected before the engine runs
        @test_throws ErrorException FUSE.ActorLocking(dd, act;
            task             = :calc_prob,
            control_type     = :EF,
            grid_size        = N_grid,
            overwrite_params = true,
            prob_method      = :conv,
            conv_window_C1   = 4)

        # :conv and :kde each build their own model type from the same grid
        for (method, T) in ((:conv, FUSE.ModeLocking.ConvProbModel),
                            (:kde,  FUSE.ModeLocking.KDEProbModel))
            actor = FUSE.ActorLocking(dd, act;
                        task             = :calc_prob,
                        control_type     = :EF,
                        grid_size        = N_grid,
                        overwrite_params = true,
                        prob_method      = method,
                        conv_window_C1   = 3,
                        conv_window_C2   = 3)
            @test actor.results.prob isa T
            C1 = actor.ode_params.Control1[1]
            C2 = actor.ode_params.Control2[1]
            @test 0.0 ≤ actor.results.prob(C1, C2) ≤ 1.0
        end
    end

end  # @testset "ActorLocking"
