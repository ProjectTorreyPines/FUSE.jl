#= ============================================================ =#
#  Test for the fuse29 pedestal predictor wrapper.                 #
#                                                                  #
#  Activates the FUSE project and exercises nn_predictor.jl in     #
#  isolation (no actors, no dd, no ZMQ). Offline-safe: the parts   #
#  that need the ONNX bundle are skipped with a warning when it is #
#  not on disk (set FUSE_PEDESTAL_NN_DIR, or let auto-fetch run).  #
#                                                                  #
#  Run with:                                                       #
#    cd FUSE && julia --project test/test_pedestal_nn_predictor.jl #
#= ============================================================ =#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using Statistics: mean, std
import NPZ

const __FUSE__ = abspath(joinpath(@__DIR__, ".."))
include(joinpath(__FUSE__, "src", "actors", "pedestal", "nn_predictor.jl"))

const REFERENCE_NPZ = joinpath(__FUSE__, "data", "pedestal_nn", "fuse29_reference_shot.npz")

@testset "fuse29 contract constants" begin
    @test length(FUSE29_ACTUATOR_NAMES) == FUSE29_N_ACTUATORS == 29
    @test length(FUSE29_ACTUATOR_UNITS) == 29
    @test length(FUSE29_HEAD_NAMES) == FUSE29_N_HEADS == 9
    @test length(FUSE29_HEAD_UNITS) == 9
    @test FUSE29_ACTUATOR_NAMES[1] == "pinj"
    @test FUSE29_ACTUATOR_NAMES[end] == "bt"
    @test FUSE29_HEAD_NAMES[1] == "ne"
    @test FUSE29_HEAD_NAMES[end] == "hmode"
    @test all(h in FUSE29_HEAD_NAMES for h in FUSE29_GATED_HEADS)
    @test !("ne" in FUSE29_GATED_HEADS) && !("hmode" in FUSE29_GATED_HEADS)
    # the old 32-channel plasma-response inputs must be gone
    @test !any(n in FUSE29_ACTUATOR_NAMES for n in ("pohm", "ip", "ipspr15v"))
    # order survives the migration from the 32-channel ResNet list (MIGRATION_FROM_RESNET.md §4)
    resnet32 = vcat(["pohm", "pinj", "tinj", "ech_total"],
                    ["f$(i)a" for i in 1:9], ["f$(i)b" for i in 1:9],
                    ["ecoila", "ecoilb"], ["gas$(c)_cal" for c in "abcde"],
                    ["ip", "ipspr15v", "bt"])
    @test filter(n -> !(n in ("pohm", "ip", "ipspr15v")), resnet32) == collect(FUSE29_ACTUATOR_NAMES)
end

@testset "fuse29 file manifest & auto-fetch surface" begin
    files = pedestal_nn_files()
    @test length(files) == 1 + length(FUSE29_BUNDLE_FILES)
    @test files[1] == "manifest.json"
    @test all(startswith(f, "s0/") for f in files[2:end])
    @test "s0/fuse29_step.onnx" in files && "s0/fuse29_seq.onnx" in files && "s0/model_config.json" in files
    @test pedestal_nn_files(; seed=3)[2] == "s3/fuse29_step.onnx"
    @test PEDESTAL_NN_HF_REPO_DEFAULT == "SCS-Lab/FUSE29-Pedestal-model"

    mktempdir() do tmp
        @test !pedestal_nn_dir_complete(tmp)
        # zero-byte stragglers count as missing
        for rel in pedestal_nn_files()
            mkpath(dirname(joinpath(tmp, rel)))
            touch(joinpath(tmp, rel))
        end
        @test !pedestal_nn_dir_complete(tmp)
        for rel in pedestal_nn_files()
            write(joinpath(tmp, rel), "x")
        end
        @test pedestal_nn_dir_complete(tmp)
        @test !pedestal_nn_dir_complete(tmp; seed=1)
        # a bare bundle dir is also "complete"
        @test pedestal_nn_dir_complete(joinpath(tmp, "s0"))
    end

    # auto-fetch opt-out never touches the network
    withenv(PEDESTAL_NN_AUTOFETCH_ENV => "0") do
        mktempdir() do tmp
            @test _maybe_autofetch!(tmp) == false
            @test !pedestal_nn_dir_complete(tmp)
        end
    end

    # env-var resolution
    withenv(PEDESTAL_NN_ENV => "/some/dir", PEDESTAL_NN_SEED_ENV => "5") do
        @test resolve_pedestal_nn_dir() == "/some/dir"
        @test _default_seed() == 5
    end
    withenv(PEDESTAL_NN_ENV => nothing, PEDESTAL_NN_SEED_ENV => nothing) do
        @test endswith(resolve_pedestal_nn_dir(), joinpath("pedestal-predictor-mamba-onnx-fuse", "artifacts"))
        @test _default_seed() == 0
    end
    @test resolve_pedestal_nn_dir(; onnx_dir="/explicit") == "/explicit"
end

# ── everything below needs the bundle on disk ────────────────────────────────
onnx_dir = resolve_pedestal_nn_dir()
have_bundle = withenv(PEDESTAL_NN_AUTOFETCH_ENV => "0") do
    pedestal_nn_dir_complete(onnx_dir)
end
if !have_bundle
    @warn "Skipping fuse29 inference tests: ONNX bundle not found at $onnx_dir (set FUSE_PEDESTAL_NN_DIR or run FUSE.download_pedestal_nn!)"
else

nn = load_pedestal_nn()

@testset "fuse29 load & contract assertions" begin
    @test load_pedestal_nn() === nn                  # process cache
    @test nn.actuator_names == collect(FUSE29_ACTUATOR_NAMES)
    @test nn.head_names == collect(FUSE29_HEAD_NAMES)
    @test nn.period_s ≈ 0.05
    @test nn.rho_ref ≈ 0.85
    @test nn.seq_len == 256
    @test nn.conv_shape == (6, 1, 1536, 3)
    @test nn.ssm_shape == (6, 1, 8, 128, 128)
    @test nn.actuator_ranges !== nothing
    @test ONNXRunTime.input_names(nn.step) == collect(FUSE29_STEP_INPUTS)
    @test ONNXRunTime.output_names(nn.step) == collect(FUSE29_STEP_OUTPUTS)

    s0 = fuse29_init_state(nn)
    @test size(s0.conv) == reverse(nn.conv_shape) && size(s0.ssm) == reverse(nn.ssm_shape)   # ORT-native layout
    @test all(iszero, s0.conv) && all(iszero, s0.ssm)
    s1 = copy(s0)
    @test s1 !== s0 && s1.conv !== s0.conv

    med = fuse29_median_actuators(nn)
    @test length(med) == 29
    @test isempty(fuse29_check_units(nn, med))
    @test isempty(fuse29_in_distribution(nn, med))
    # a factor-of-1e6 unit error on pinj must be caught
    bad = copy(med); bad[1] = 3.0e6
    w = fuse29_check_units(nn, bad)
    @test length(w) == 1 && startswith(w[1], "pinj")
    @test length(fuse29_in_distribution(nn, bad)) == 1
end

@testset "fuse29 reference shot replay (both graphs)" begin
    # actuator_names/head_names/readme are numpy '<U' unicode arrays, which NPZ.jl
    # cannot parse; read only the numeric keys.
    ref = NPZ.npzread(REFERENCE_NPZ, ["actuators", "predictions", "seed", "period_ms", "shot"])
    A = Float32.(ref["actuators"])
    P = Float32.(ref["predictions"])
    @test size(A, 2) == 29 && size(P, 2) == 9 && size(A, 1) == size(P, 1)
    @test Float64(ref["period_ms"]) ≈ nn.period_s * 1000
    @test isempty(fuse29_check_units(nn, A))

    K = size(A, 1)
    R = fuse29_rollout(nn, A)
    S = fuse29_sequence(nn, A)
    @test size(R) == (K, 9) == size(S)
    @test maximum(abs, R .- S) < 1e-3                      # step graph == sequence graph
    if Int(ref["seed"]) == nn.seed
        @test maximum(abs, R .- P) < 1e-3                  # matches the shipped reference (seed 0)
        @test maximum(abs, S .- P) < 1e-3
    else
        @warn "bundle is seed $(nn.seed), reference is seed $(ref["seed"]); skipping exact comparison"
    end
    # right-padding is exact: predictions for the real ticks do not change
    S2 = fuse29_sequence(nn, A[1:K-10, :])
    @test maximum(abs, S2 .- S[1:K-10, :]) < 1e-5
    @test_throws ErrorException fuse29_sequence(nn, zeros(Float32, nn.seq_len + 1, 29))

    hm = R[:, end]
    @test all(0 .<= hm .<= 1)
    @test any(hm .< 0.5) && any(hm .> 0.5)                 # the shot has an L-H transition
    @test 0.5 < mean(R[:, 1]) < 10                          # ne, 1e19 m^-3
    @test 0.05 < mean(R[:, 2]) < 3                          # te_ped, keV
end

@testset "fuse29_step is functional; commit semantics" begin
    ref = NPZ.npzread(REFERENCE_NPZ, ["actuators"])
    A = Float32.(ref["actuators"])
    s = fuse29_init_state(nn)
    for k in 1:20
        _, s = fuse29_step(nn, A[k, :], s)
    end
    conv_before = copy(s.conv); ssm_before = copy(s.ssm)
    p1, s1 = fuse29_step(nn, A[21, :], s)
    @test s.conv == conv_before && s.ssm == ssm_before     # input state untouched
    # re-stepping N times from the same committed state == one tick
    for _ in 1:4
        p, sn = fuse29_step(nn, A[21, :], s)
        @test p == p1
        @test sn.ssm == s1.ssm
    end
    @test_throws ErrorException fuse29_step(nn, A[21, 1:28], s)
    @test_throws ErrorException fuse29_step(nn, [NaN32; A[21, 2:end]], s)

    nt = fuse29_named(p1)
    @test keys(nt) == Symbol.(FUSE29_HEAD_NAMES)
    pr = fuse29_prediction(nn, p1, false; time=1.0, n_ticks=1)
    @test isnan(pr.neped_prmtan) && isnan(pr.rho_sym) && isfinite(pr.ne) && isfinite(pr.te_ped)
    @test pr.hmode_prob == Float64(p1[end]) && pr.is_h_mode == false && pr.rho_ref == nn.rho_ref
    pr = fuse29_prediction(nn, p1, true)
    @test isfinite(pr.neped_prmtan) && pr.is_h_mode
end

@testset "Fuse29Tracker: one tick per 50 ms, hold, catch-up, re-run, rewind" begin
    med = fuse29_median_actuators(nn)
    calls = Float64[]
    act = t -> (push!(calls, t); med)
    period = nn.period_s

    tr = Fuse29Tracker(nn, 1.0 - period)           # first call at t=1.0 -> exactly one tick
    p1, h1, n1 = fuse29_advance!(tr, nn, 1.0, act)
    @test n1 == 1 && calls ≈ [1.0]
    @test tr.pending_time ≈ 1.0 && tr.committed_time ≈ 1.0 - period

    # re-run of the same step: nothing committed, re-stepped from committed
    empty!(calls)
    p1b, _, n1b = fuse29_advance!(tr, nn, 1.0, act)
    @test n1b == 1 && p1b == p1 && tr.committed_time ≈ 1.0 - period && tr.n_total == 2

    # host steps faster than 50 ms: zero-order hold, no model tick
    empty!(calls)
    p2, _, n2 = fuse29_advance!(tr, nn, 1.02, act)
    @test n2 == 0 && isempty(calls) && p2 == p1
    @test tr.committed_time ≈ 1.0                  # the 1.0 tick got committed when time advanced

    # host steps slower than 50 ms: catch-up with one call per 50 ms sub-interval
    empty!(calls)
    p3, _, n3 = fuse29_advance!(tr, nn, 1.10, act)
    @test n3 == 2 && calls ≈ [1.05, 1.10]
    @test tr.pending_time ≈ 1.10

    # consistency: the tracker's trajectory equals a plain rollout on the same grid
    R = fuse29_rollout(nn, repeat(med', 3, 1))
    @test maximum(abs, p3 .- R[3, :]) < 1e-5

    # rewind is an error at this level (the actor resets the tracker)
    @test_throws ErrorException fuse29_advance!(tr, nn, 0.5, act)
    @test_throws ErrorException fuse29_advance!(tr, nn, 1.2, act; hmode_enter=0.3, hmode_exit=0.7)

    # Schmitt trigger: enter above `enter`, stay until below `exit`
    tr2 = Fuse29Tracker(nn, 0.0)
    hs = Bool[]
    for k in 1:40
        _, h, _ = fuse29_advance!(tr2, nn, k * period, act; hmode_enter=0.7, hmode_exit=0.3)
        push!(hs, h)
    end
    @test hs[1] == false                            # first tick from zero state is L-mode-ish
    @test any(hs)                                   # median actuators reach H-mode within 2 s
    @test count(i -> hs[i] != hs[i-1], 2:length(hs)) <= 2
end

end # have_bundle
