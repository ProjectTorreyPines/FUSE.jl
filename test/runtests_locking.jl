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

    # STRUCTURE: one testset group per control_type. :solve_system is the expensive
    # step (grid_size^2 ODE trajectories), and :calc_prob / :eval_prob / :transfer_learning
    # only CONSUME its output — they load from disk, keyed by control_type
    # (save_ode_results / load_ode_results!). So each group solves ONCE at group scope
    # and every testset below it reuses that solve. Ordering inside a group matters:
    # calc_prob writes nn_model_<control_type>.bson, which eval_prob and
    # transfer_learning both read.
    #
    # The suite runs with overwrite_params=true throughout, so the geometry is pinned
    # to PoP2024 (r0 = 1 m, r_w = 1.0, r_c = 1.25) and so are mu and Inertia. That
    # keeps the tests off dd.core_sources carrying a torque, which the actor errors on
    # rather than falling back. The geometry testset at the bottom covers the
    # dd-derived path (overwrite_params=false) explicitly, and cheaply, via task=:bounds.

    # Small grid and minimal NN for test speed
    N_grid  = 10
    fast_nn = FUSE.NNparams(hidden_sizes=[10], n_epochs=10, batch_size=8)

    # Control2_min/max have no default — the field is Gauss for :EF, Δ_RW for :LinStab
    # and α for :NLsaturation, so the actor refuses to guess. Each group names its own.
    ef_range = (Control2_min=0.01, Control2_max=10.0)   # Gauss

    # ═══ control_type = :EF — error field ════════════════════════════════════
    @testset verbose = true "EF control" begin
        # The one EF solve. Every testset below reads it back from disk.
        actor_EF = FUSE.ActorLocking(dd, act;
                       task             = :solve_system,
                       control_type     = :EF,
                       grid_size        = N_grid,
                       overwrite_params = true,
                       ef_range...)
        saved_ode = copy(actor_EF.results.ode_sols)

        @testset "single_case, one ODE trajectory" begin
            # Solve one ODE trajectory only — no grid, no NN, no stored results
            actor = FUSE.ActorLocking(dd, act;
                        task             = :single_case,
                        op_C2            = 5.0,
                        overwrite_params = true,
                        ef_range...)
            @test actor.results === nothing
        end

        @testset "solve_system, grid and bifurcation bounds" begin
            r = actor_EF.results

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

        @testset "calc_prob, retrain NN without re-solving" begin
            # Fresh actor: results === nothing, so it loads the grid the group
            # solve checkpointed and retrains on it
            actor_prob = FUSE.ActorLocking(dd, act;
                             task             = :calc_prob,
                             control_type     = :EF,
                             grid_size        = N_grid,
                             overwrite_params = true,
                             ef_range...,
                             nn_params        = fast_nn)
            @test actor_prob.results !== nothing
            @test actor_prob.results.prob !== nothing
            # ODE grid must be identical to what was saved — not re-computed
            @test actor_prob.results.ode_sols ≈ saved_ode
        end

        @testset "calc_prob, prob_method engines" begin
            # an even window has no centre point — rejected before the engine runs
            @test_throws ErrorException FUSE.ActorLocking(dd, act;
                task             = :calc_prob,
                control_type     = :EF,
                grid_size        = N_grid,
                overwrite_params = true,
                ef_range...,
                prob_method      = :conv,
                conv_window_C1   = 4)

            # :conv and :kde each build their own model type from the same grid.
            # They save under their own filenames, so nn_model_EF.bson survives for
            # the eval_prob and transfer_learning testsets below.
            for (method, T) in ((:conv, FUSE.ModeLocking.ConvProbModel),
                                (:kde,  FUSE.ModeLocking.KDEProbModel))
                actor = FUSE.ActorLocking(dd, act;
                            task             = :calc_prob,
                            control_type     = :EF,
                            grid_size        = N_grid,
                            overwrite_params = true,
                            ef_range...,
                            prob_method      = method,
                            conv_window_C1   = 3,
                            conv_window_C2   = 3)
                @test actor.results.prob isa T
                C1 = actor.ode_params.Control1[1]
                C2 = actor.ode_params.Control2[1]
                @test 0.0 ≤ actor.results.prob(C1, C2) ≤ 1.0
            end
        end

        @testset "eval_prob, probability at chosen times" begin
            # op_times must be ≥ dd.global_time: the group solve already appended an
            # mhd_linear slice there, and resize! rejects a time below the last one
            t_ref    = dd.global_time
            op_times = t_ref .+ [0.0, 1e-3, 2e-3]

            # A scalar op_C2 broadcasts across every op_time
            actor = FUSE.ActorLocking(dd, act;
                        task             = :eval_prob,
                        control_type     = :EF,
                        grid_size        = N_grid,
                        overwrite_params = true,
                        ef_range...,
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

            # The actor asserts NO limits. It reports what it measured — stability_metric
            # on mhd_linear — and a limit model in limit_models.jl owns the thresholds and
            # writes limits.model, the same split ActorVerticalStability uses.
            @test isempty(dd.limits.model)

            # A vector op_C2 must match op_times one-for-one. overwrite_params=true so
            # the ErrorException can only be the length mismatch, not the geometry.
            @test_throws ErrorException FUSE.ActorLocking(dd, act;
                task             = :eval_prob,
                control_type     = :EF,
                grid_size        = N_grid,
                overwrite_params = true,
                ef_range...,
                op_times         = op_times,
                op_C2            = [1.0, 5.0])
        end

        @testset "transfer_learning, fine-tune on focused grid" begin
            # The base model is nn_model_EF.bson, written by the calc_prob testset above.
            # Fine-tune (last layer only) on a focused/sparser control-space sweep.
            tl_nn = FUSE.NNparams(hidden_sizes=[10], n_epochs=5, batch_size=8)
            actor_tl = FUSE.ActorLocking(dd, act;
                           task             = :transfer_learning,
                           control_type     = :EF,
                           grid_size        = 5,
                           overwrite_params = true,
                           ef_range...,
                           nn_params        = tl_nn,
                           Control1_min     = 1.0,
                           Control1_max     = 3.0)
            r = actor_tl.results

            @test r !== nothing
            @test size(r.ode_sols, 1) == 5^2
            @test r.prob isa FUSE.LockingNNModel
            C1 = actor_tl.ode_params.Control1[1]
            C2 = actor_tl.ode_params.Control2[1]
            @test 0.0 ≤ r.prob(C1, C2) ≤ 1.0
        end
    end

    # ═══ control_type = :LinStab — stability index ═══════════════════════════
    @testset verbose = true "LinStab control" begin
        # Control2_min/max are Δ_RW for :LinStab and must be < 0; the defaults
        # are positive, so they have to be given explicitly — to every task, since
        # the control grid is rebuilt before the task branch runs
        drw_range = (Control2_min=-3.5, Control2_max=-0.05)

        # The one LinStab solve, plus the model the eval testset evaluates
        actor_LS = FUSE.ActorLocking(dd, act;
                       task             = :solve_system,
                       control_type     = :LinStab, error_field=10.0,
                       grid_size        = N_grid,
                       overwrite_params = true,
                       drw_range...)
        FUSE.ActorLocking(dd, act;
            task             = :calc_prob,
            control_type     = :LinStab, error_field=10.0,
            grid_size        = N_grid,
            overwrite_params = true,
            drw_range...,
            nn_params        = fast_nn)

        # task=:bounds solves no ODEs — one actor serves both geometry testsets
        actor_bounds = FUSE.ActorLocking(dd, act;
                           task             = :bounds,
                           control_type     = :LinStab, error_field=10.0,
                           grid_size        = N_grid,
                           overwrite_params = true,
                           drw_range...)

        @testset "solve_system, grid over Δ_RW" begin
            r = actor_LS.results
            @test r !== nothing
            @test size(r.ode_sols, 1) == N_grid^2
            # NN training is deferred to task=:calc_prob — :solve_system leaves prob unset
            @test r.prob === nothing
        end

        @testset "eval_prob, br to Δt inversion" begin
            # later than the times the :EF eval_prob testset already appended
            t_ref    = dd.global_time
            op_times = t_ref .+ [3e-3, 4e-3]

            # op_C2 is the n=1 br amplitude in Gauss here, one per time
            actor = FUSE.ActorLocking(dd, act;
                        task             = :eval_prob,
                        control_type     = :LinStab, error_field=10.0,
                        grid_size        = N_grid,
                        overwrite_params = true,
                        drw_range...,
                        op_times         = op_times,
                        op_C2            = [10.0, 20.0])

            ev = actor.eval
            @test ev !== nothing
            @test ev.C2 == [10.0, 20.0]            # stored in user units, not Δt

            op = actor.ode_params
            # error_field is derived from par.error_field (Gauss), never accumulated
            @test op.error_field ≈ actor.par.error_field * 1e-4 / actor.par.mag_perturbation_amplitude *
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

        @testset "bounds, hysteresis onset without solving the ODEs" begin
            # calculate_bifurcation_bounds reads only geometry and the control grid
            b = actor_bounds.bounds
            @test b !== nothing
            @test actor_bounds.results === nothing     # nothing was solved
            @test length(b.C1) == N_grid
            @test size(b.map) == (N_grid, N_grid)
            @test length(b.C2_onset) == length(b.C1) == length(b.C2_onset_user)
            # NL_saturation off ⇒ α=0 ⇒ smooth discriminant ⇒ onset interpolated
            # between grid columns rather than bracketed by them
            @test b.bracketed == false

            # at least some rotations must resolve an onset, else the rest is vacuous
            keep = findall(isfinite, b.C2_onset_user)
            @test !isempty(keep)

            # C1_user is kHz; the onset in user units is a positive amplitude
            @test all(>(0), b.C2_onset_user[keep])
            @test issorted(b.C1_user)

            # Monotone increasing in rotation — faster rotation tolerates a larger island.
            # The inverse lookup (min rotation for a given amplitude) depends on this.
            @test issorted(b.C2_onset_user[keep])

            # :bounds asserts no limits either
            @test isempty(dd.limits.model)
        end

        @testset "bounds, grid-bracketed onset under NL saturation" begin
            # NL_saturation=true leaves ode_params.saturation_param at its nonzero
            # default instead of zeroing it, so calculate_bifurcation_bounds takes the
            # α≠0 path: per-point root counting returning ±1 instead of the cubic
            # discriminant. :LinStab is the control_type where this flag actually
            # selects the path — for :NLsaturation α is the swept axis regardless.
            actor_nl = FUSE.ActorLocking(dd, act;
                           task             = :bounds,
                           control_type     = :LinStab, error_field=10.0,
                           grid_size        = N_grid,
                           overwrite_params = true,
                           NL_saturation    = true,
                           drw_range...)

            b = actor_nl.bounds
            @test b !== nothing
            @test actor_nl.results === nothing
            @test actor_nl.ode_params.saturation_param != 0.0   # not zeroed ⇒ α≠0 path
            @test b.bracketed == true

            # bb carries the ±1 flag, not a discriminant: no other values may appear
            @test all(v -> v == 1.0 || v == -1.0, b.map)

            keep = findall(isfinite, b.C2_onset)
            @test !isempty(keep)

            # A ±1 flag has no within-cell information, so the onset cannot be
            # interpolated — every resolved value must land exactly on a scanned column
            c2_axis = unique(actor_nl.ode_params.Control2)
            @test all(b.C2_onset[k] ∈ c2_axis for k in keep)

            # Monotone in rotation still, but only weakly: neighbouring rotations can
            # bracket into the same column, so ties are expected here and issorted
            # is the right (non-strict) check
            @test issorted(b.C2_onset_user[findall(isfinite, b.C2_onset_user)])
        end

        @testset "br_drw, PoP2024 quadratic and its alpha=0 limit" begin
            op = actor_bounds.ode_params

            # α = 0 must reproduce the closed form ψt = P_rw/Δ_RW exactly
            P_rw = op.l21 * op.l32 * FUSE._eps_ref(actor_bounds) / op.DeltaW
            b0   = actor_bounds.par.mag_perturbation_amplitude
            for drw in (-2.0, -0.5)
                br = FUSE.br_drw(actor_bounds; direction=:forward, drw, alpha=0.0, verbose=false)
                ψ  = (br / (b0 * 1e4)) * op.rat_surface / actor_bounds.par.m_pol
                @test ψ ≈ P_rw / drw
            end

            # α ≠ 0 takes the upper root of the quadratic, and forward/backward invert
            for α in (0.0, 0.05, 0.2, 0.5), drw in (-2.0, -0.5)
                br   = FUSE.br_drw(actor_bounds; direction=:forward, drw, alpha=α, verbose=false)
                back = FUSE.br_drw(actor_bounds; direction=:backward, br_Gauss=br, alpha=α, verbose=false)
                @test back ≈ drw rtol=1e-10
            end

            # saturation reduces the amplitude: ψ(α) is strictly decreasing
            brs = [FUSE.br_drw(actor_bounds; direction=:forward, drw=-2.0, alpha=α, verbose=false)
                   for α in (0.0, 0.1, 0.3, 0.6)]
            @test issorted(brs; rev=true)
        end
    end

    # ═══ control_type = :NLsaturation ════════════════════════════════════════
    @testset "NLsaturation control, solve_system over the alpha grid" begin
        actor = FUSE.ActorLocking(dd, act;
                    task             = :solve_system,
                    control_type     = :NLsaturation, error_field=10.0,
                    NL_saturation    = true,
                    grid_size        = N_grid,
                    overwrite_params = true,
                    Control2_min     = 0.05,
                    Control2_max     = 0.5)
        r = actor.results
        @test r !== nothing
    end

    # ═══ Control2 bounds have no default ═════════════════════════════════════
    @testset "Control2_min/max are required, not defaulted" begin
        # Control2 is Gauss / Δ_RW / α depending on control_type, so a default that
        # suits one silently mis-scans the other two. Every control_type must refuse.
        for ct in (:EF, :LinStab, :NLsaturation)
            @test_throws ErrorException FUSE.ActorLocking(dd, act;
                task=:bounds, control_type=ct, grid_size=N_grid, overwrite_params=true,
                error_field=10.0)
        end

        # error_field has no default either, but only :LinStab/:NLsaturation read it
        @test_throws ErrorException FUSE.ActorLocking(dd, act;
            task=:bounds, control_type=:LinStab, grid_size=N_grid,
            overwrite_params=true, Control2_min=-3.5, Control2_max=-0.05)

        # one bound alone is not enough either
        @test_throws ErrorException FUSE.ActorLocking(dd, act;
            task=:bounds, control_type=:LinStab, error_field=10.0, grid_size=N_grid,
            overwrite_params=true, Control2_min=-3.5)

        # and the :LinStab sign check still fires once both are given, so the NaN
        # guard has not displaced the validation that follows it
        @test_throws ErrorException FUSE.ActorLocking(dd, act;
            task=:bounds, control_type=:LinStab, error_field=10.0, grid_size=N_grid,
            overwrite_params=true, Control2_min=0.01, Control2_max=10.0)
    end

    # ═══ geometry — r0 comes from dd unless overwrite_params pins it ═════════
    @testset "geometry, r0 from dd.equilibrium and the pinned radii" begin
        a_minor = dd.equilibrium.time_slice[].boundary.minor_radius
        @test a_minor > 0

        # default (NaN) resolves to the minor radius; an explicit value wins;
        # overwrite_params pins r0=1.0 for PoP2024 regardless of either — and, since
        # the radii are in metres, also pins r_w=1.0 and r_c=1.25.
        # task=:bounds so this exercises the dd-derived geometry without an ODE solve.
        actor_dd = FUSE.ActorLocking(dd, act;
                       task=:bounds, control_type=:LinStab, error_field=10.0, grid_size=N_grid,
                       overwrite_params=false, Control2_min=-3.5, Control2_max=-0.05)
        @test FUSE._length_scale(dd, actor_dd.par) ≈ a_minor
        # radii are metres: the derived dimensionless value tracks r0
        @test actor_dd.ode_params.control_surf ≈ actor_dd.par.control_surf_radius / a_minor

        actor_fixed = FUSE.ActorLocking(dd, act;
                          task=:bounds, control_type=:LinStab, error_field=10.0, grid_size=N_grid,
                          overwrite_params=false, r0=0.75,
                          Control2_min=-3.5, Control2_max=-0.05)
        @test FUSE._length_scale(dd, actor_fixed.par) == 0.75

        actor_ow = FUSE.ActorLocking(dd, act;
                       task=:bounds, control_type=:LinStab, error_field=10.0, grid_size=N_grid,
                       overwrite_params=true, r0=0.75,
                       res_wall_radius=0.6, control_surf_radius=0.9,
                       Control2_min=-3.5, Control2_max=-0.05)
        @test FUSE._length_scale(dd, actor_ow.par) == 1.0
        # overwrite_params pins the whole geometry, not just r0: supplied radii ignored
        @test actor_ow.ode_params.control_surf == 1.25
    end

end  # @testset "ActorLocking"
