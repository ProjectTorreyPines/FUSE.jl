#= ============================================================ =#
#  Test for the fuse_efit pedestal predictor wrapper.              #
#                                                                  #
#  Activates the FUSE project and exercises nn_predictor.jl +      #
#  nn_predictor_efit.jl in isolation (no actors, no dd, no ZMQ).   #
#  Offline-safe: the parts that need the ONNX bundle are skipped   #
#  with a warning when it is not on disk (set                      #
#  FUSE_PEDESTAL_EFIT_NN_DIR, or let auto-fetch run).              #
#                                                                  #
#  Run with:                                                       #
#    cd FUSE && julia --project test/test_pedestal_nn_efit.jl      #
#= ============================================================ =#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
import NPZ

const __FUSE__ = abspath(joinpath(@__DIR__, ".."))
include(joinpath(__FUSE__, "src", "actors", "pedestal", "nn_predictor.jl"))
include(joinpath(__FUSE__, "src", "actors", "pedestal", "nn_predictor_efit.jl"))

const REFERENCE_NPZ = joinpath(__FUSE__, "data", "pedestal_nn", "fuse_efit_reference_shot.npz")

@testset "fuse_efit contract constants" begin
    @test length(FUSE_EFIT_INPUT_NAMES) == FUSE_EFIT_N_INPUTS == 25
    @test length(FUSE_EFIT_INPUT_UNITS) == 25
    @test length(FUSE_EFIT_INDICATOR_NAMES) == FUSE_EFIT_N_INDICATORS == 5
    @test length(FUSE_EFIT_HEAD_NAMES) == FUSE_EFIT_N_HEADS == 9
    @test length(FUSE_EFIT_HEAD_UNITS) == 9
    @test FUSE_EFIT_INPUT_NAMES[1] == "pinj" && FUSE_EFIT_INPUT_NAMES[10] == "ip" && FUSE_EFIT_INPUT_NAMES[end] == "vloop"
    # the 14 shape channels are columns 11..24 (0-based 10..23), the block `efit_valid` governs
    @test collect(FUSE_EFIT_INPUT_NAMES[11:24]) == collect(FUSE_EFIT_SHAPE_NAMES)
    # the coil currents of fuse29 are gone; the actuators the two models share kept their names
    @test !any(startswith(n, "f") && length(n) == 3 for n in FUSE_EFIT_INPUT_NAMES)
    @test !("ecoila" in FUSE_EFIT_INPUT_NAMES) && !("ech_total" in FUSE_EFIT_INPUT_NAMES)
    @test all(n in FUSE29_ACTUATOR_NAMES for n in ("pinj", "tinj", "gasa_cal", "gase_cal", "bt"))
    @test all(n in FUSE_EFIT_INPUT_NAMES for n in ("pinj", "tinj", "gasa_cal", "gase_cal", "bt"))
    # same head layout as fuse29: the continuous heads and hmode sit in the same columns
    @test FUSE_EFIT_HEAD_NAMES[1:4] == FUSE29_HEAD_NAMES[1:4] && FUSE_EFIT_HEAD_NAMES[end] == FUSE29_HEAD_NAMES[end] == "hmode"
    @test all(h in FUSE_EFIT_HEAD_NAMES for h in FUSE_EFIT_GATED_HEADS)
    @test !("ne" in FUSE_EFIT_GATED_HEADS) && !("hmode" in FUSE_EFIT_GATED_HEADS)
end

@testset "fuse_efit file manifest & auto-fetch surface" begin
    files = pedestal_efit_nn_files()
    @test files[1:2] == ["manifest.json", "SHA256SUMS"]
    @test length(files) == 2 + 6
    @test "era/single/fuse_efit_single_step.onnx" in files && "era/single/fuse_efit_single_sequence.onnx" in files
    @test "era/averaged/fuse_efit_averaged_step.onnx" in pedestal_efit_nn_files(; form=:averaged)
    @test_throws ErrorException fuse_efit_bundle_files(:nope)
    @test PEDESTAL_EFIT_NN_HF_REPO_DEFAULT == "SCS-Lab/FUSE-PED"

    mktempdir() do tmp
        @test _fuse_efit_find_bundle_dir(tmp, :single) === nothing
        for rel in pedestal_efit_nn_files()
            mkpath(dirname(joinpath(tmp, rel)))
            touch(joinpath(tmp, rel))               # zero-byte stragglers count as missing
        end
        @test _fuse_efit_find_bundle_dir(tmp, :single) === nothing
        for rel in pedestal_efit_nn_files()
            write(joinpath(tmp, rel), "x")
        end
        # the download root, the `era` folder and the bare bundle folder are all accepted
        @test _fuse_efit_find_bundle_dir(tmp, :single) == joinpath(tmp, "era", "single")
        @test _fuse_efit_find_bundle_dir(joinpath(tmp, "era"), :single) == joinpath(tmp, "era", "single")
        @test _fuse_efit_find_bundle_dir(joinpath(tmp, "era", "single"), :single) == joinpath(tmp, "era", "single")
        @test _fuse_efit_find_bundle_dir(tmp, :averaged) === nothing
        # size check against a manifest
        write(joinpath(tmp, "manifest.json"), """{"files": {"era/single/norms.json": {"bytes": 1}, "era/single/model_config.json": {"bytes": 7}}}""")
        msgs = fuse_efit_verify_sizes(tmp)
        @test length(msgs) == 1 && occursin("model_config.json", msgs[1])
    end

    # auto-fetch opt-out never touches the network
    withenv(PEDESTAL_NN_AUTOFETCH_ENV => "0") do
        mktempdir() do tmp
            @test _fuse_efit_maybe_autofetch!(tmp, :single) == false
            @test_throws ErrorException load_pedestal_efit_nn(; onnx_dir=tmp)
        end
    end
    withenv(PEDESTAL_EFIT_NN_ENV => "/some/where") do
        @test resolve_pedestal_efit_nn_dir() == "/some/where"
        @test resolve_pedestal_efit_nn_dir(; onnx_dir="/else") == "/else"
    end
end

have_bundle = withenv(PEDESTAL_NN_AUTOFETCH_ENV => "0") do
    _fuse_efit_find_bundle_dir(abspath(resolve_pedestal_efit_nn_dir()), :single) !== nothing
end
if !have_bundle
    @warn "fuse_efit bundle not on disk ($(resolve_pedestal_efit_nn_dir())); skipping the model tests of test_pedestal_nn_efit"
else

nn = load_pedestal_efit_nn()

@testset "fuse_efit load: contract asserted against the bundle" begin
    @test nn.form == :single
    @test nn.input_names == collect(FUSE_EFIT_INPUT_NAMES)
    @test nn.indicator_names == collect(FUSE_EFIT_INDICATOR_NAMES)
    @test nn.head_names == collect(FUSE_EFIT_HEAD_NAMES)
    @test nn.period_s == 0.05 && nn.rho_ref == 0.85 && nn.seq_len == 256
    @test nn.conv_shape == (6, 1, 1792, 3) && nn.ssm_shape == (6, 1, 12, 128, 64)
    @test nn.hmode_balanced == 0.54 && nn.hmode_conservative == 0.50
    @test load_pedestal_efit_nn() === nn                       # cached
    @test isempty(fuse_efit_verify_sizes(abspath(resolve_pedestal_efit_nn_dir())))
    s0 = fuse_efit_init_state(nn)
    @test size(s0.conv) == reverse(nn.conv_shape) && size(s0.ssm) == reverse(nn.ssm_shape)
    @test all(iszero, s0.conv) && all(iszero, s0.ssm)
end

ref = NPZ.npzread(REFERENCE_NPZ, ["inputs", "indicators", "available", "expected_single", "tolerance", "period_ms", "shot"])
X, I, av, E, tol = ref["inputs"], ref["indicators"], ref["available"], ref["expected_single"], ref["tolerance"]
K = size(X, 1)

@testset "fuse_efit reference shot $(Int(ref["shot"])): step and sequence graphs reproduce the model authors' outputs" begin
    @test size(X) == (K, 25) && size(I) == (K, 5) && size(E) == (K, 9) && length(tol) == 9
    @test Float64(ref["period_ms"]) / 1000 == nn.period_s
    R = fuse_efit_rollout(nn, X, I, av)
    S = fuse_efit_sequence(nn, X, I, av)
    for j in 1:9
        @test maximum(abs.(R[:, j] .- E[:, j])) <= tol[j]      # routing, units, column order, state layout
        @test maximum(abs.(S[:, j] .- E[:, j])) <= tol[j]
        @test maximum(abs.(R[:, j] .- S[:, j])) <= tol[j]      # step graph == sequence graph
    end
    @test all(0 .<= R[:, end] .<= 1)                            # hmode is a probability
    @test_throws ErrorException fuse_efit_sequence(nn, zeros(Float32, nn.seq_len + 1, 25), ones(Float32, nn.seq_len + 1, 5))
    @test_throws ErrorException fuse_efit_rollout(nn, X[:, 1:24], I, av)
end

@testset "fuse_efit: masked values cannot reach the model" begin
    # this shot has no upper strike points on most ticks: rvsou / rvsiu are NaN there with their flags at 0
    k = findfirst(i -> I[i, 4] == 0 && I[i, 5] == 0 && I[i, 1] == 1, 1:K)
    @test k !== nothing
    s0 = fuse_efit_init_state(nn)
    p_ref, _ = fuse_efit_step(nn, X[k, :], I[k, :], av, s0)
    jou = findfirst(==("rvsou"), nn.input_names)
    jiu = findfirst(==("rvsiu"), nn.input_names)
    for junk in (0.0f0, -0.89f0, 123.0f0, NaN32)
        x = copy(X[k, :])
        x[jou] = junk
        x[jiu] = junk
        p, _ = fuse_efit_step(nn, x, I[k, :], av, s0)
        @test p == p_ref
    end
    # efit_valid = 0 masks the whole shape block; vloop_valid = 0 masks the loop voltage
    ind = copy(I[k, :])
    ind[1] = 0
    ind[2] = 0
    x = copy(X[k, :])
    pa, _ = fuse_efit_step(nn, x, ind, av, s0)
    x[11:25] .= NaN32
    pb, _ = fuse_efit_step(nn, x, ind, av, s0)
    @test pa == pb && pa != p_ref
    # an unavailable channel is ignored too, and the flag itself is an input (so the output moves)
    av0 = copy(av)
    jg = findfirst(==("gase_cal"), nn.input_names)
    av0[jg] = 0
    x = copy(X[k, :])
    pc, _ = fuse_efit_step(nn, x, I[k, :], av0, s0)
    x[jg] = 1.0f6
    pd, _ = fuse_efit_step(nn, x, I[k, :], av0, s0)
    @test pc == pd && pc != p_ref
end

@testset "fuse_efit step is functional; prediction record" begin
    s0 = fuse_efit_init_state(nn)
    p1, s1 = fuse_efit_step(nn, X[1, :], I[1, :], av, s0)
    @test all(iszero, s0.conv) && all(iszero, s0.ssm)           # input state untouched
    p1b, _ = fuse_efit_step(nn, X[1, :], I[1, :], av, s0)
    @test p1 == p1b
    p2, _ = fuse_efit_step(nn, X[2, :], I[2, :], av, s1)
    @test maximum(abs.(p2 .- E[2, :]) ./ tol) <= 1
    @test_throws ErrorException fuse_efit_step(nn, X[1, 1:24], I[1, :], av, s0)
    @test_throws ErrorException fuse_efit_step(nn, X[1, :], I[1, 1:4], av, s0)

    @test keys(fuse_efit_named(p1)) == Symbol.(FUSE_EFIT_HEAD_NAMES)
    pr = fuse_efit_prediction(nn, p1, false; time=1.0, n_ticks=1)
    @test isnan(pr.ne_knee) && isnan(pr.te_knee) && isnan(pr.rho_sym) && isnan(pr.ne_top_loc)
    @test isfinite(pr.ne) && isfinite(pr.te_ped) && isfinite(pr.ti_ped) && isfinite(pr.t_rot_ped)
    @test pr.hmode_prob == Float64(p1[end]) && pr.is_h_mode == false && pr.rho_ref == 0.85 && pr.model == :fuse_efit
    pr = fuse_efit_prediction(nn, p1, true)
    @test isfinite(pr.ne_knee) && pr.is_h_mode
    # every field ActorPedestal reads from a fuse29 prediction is there
    @test all(hasproperty(pr, f) for f in (:ne, :te_ped, :ti_ped, :t_rot_ped, :hmode_prob, :is_h_mode, :rho_ref, :time, :n_ticks, :raw))
end

@testset "fuse_efit units smoke test" begin
    k = 60
    @test isempty(fuse_efit_check_units(nn, X[k, :], I[k, :], av))
    x = copy(X[k, :])
    jip = findfirst(==("ip"), nn.input_names)
    x[jip] *= 1e6                                              # amperes for mega-amperes
    w = fuse_efit_check_units(nn, x, I[k, :], av)
    @test length(w) == 1 && startswith(w[1], "ip")
    x = copy(X[k, :])
    ja = findfirst(==("aminor"), nn.input_names)
    x[ja] *= 100                                               # centimetres for metres
    @test any(startswith(m, "aminor") for m in fuse_efit_check_units(nn, x, I[k, :], av))
    ind = copy(I[k, :])
    ind[1] = 0                                                 # ... unless the block is masked
    @test !any(startswith(m, "aminor") for m in fuse_efit_check_units(nn, x, ind, av))
    @test all(e.channel in FUSE_EFIT_INPUT_NAMES for e in fuse_efit_in_distribution(nn, X[k, :], I[k, :], av))
end

@testset "PedestalNNTracker with fuse_efit: one tick per 50 ms, hold, catch-up, re-run" begin
    period = nn.period_s
    calls = Float64[]
    # tick index of time t on the reference shot's grid, starting at tick 1 for t = period
    step_at = (t, s) -> begin
        push!(calls, t)
        k = round(Int, t / period)
        fuse_efit_step(nn, X[k, :], I[k, :], av, s)
    end
    tr = PedestalNNTracker(nn, 0.0)
    @test length(tr.last) == 9
    p1, _, n1 = pedestal_nn_advance!(tr, period, period, step_at)
    @test n1 == 1 && calls ≈ [period]
    p1b, _, n1b = pedestal_nn_advance!(tr, period, period, step_at)             # re-run: re-stepped from committed
    @test n1b == 1 && p1b == p1 && tr.committed_time ≈ 0.0
    empty!(calls)
    p2, _, n2 = pedestal_nn_advance!(tr, period, period + 0.02, step_at)        # faster than 50 ms: hold
    @test n2 == 0 && isempty(calls) && p2 == p1
    p3, _, n3 = pedestal_nn_advance!(tr, period, 4 * period, step_at)           # slower: catch-up
    @test n3 == 3 && calls ≈ [2, 3, 4] .* period
    @test maximum(abs.(p3 .- E[4, :]) ./ tol) <= 1                              # same trajectory as the reference
    @test_throws ErrorException pedestal_nn_advance!(tr, period, 0.0, step_at)  # rewind is the actor's job

    # the gate with the form's own thresholds follows the reference shot's hmode probability
    tr2 = PedestalNNTracker(nn, 0.0)
    hs = Bool[]
    for k in 1:K
        _, h, _ = pedestal_nn_advance!(tr2, period, k * period, step_at; hmode_enter=nn.hmode_balanced, hmode_exit=nn.hmode_conservative)
        push!(hs, h)
    end
    @test hs[1] == (E[1, end] > nn.hmode_balanced)
    @test all(hs[k] for k in 1:K if E[k, end] > nn.hmode_balanced + 0.01)
    @test !any(hs[k] for k in 1:K if E[k, end] < nn.hmode_conservative - 0.01)
end

end # have_bundle
