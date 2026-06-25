#= ============================================================ =#
#  Smoke test for the PedestalNN ensemble wrapper.                #
#                                                                  #
#  Activates the FUSE project and exercises the nn_predictor.jl    #
#  module in isolation (no actors, no dd, no ZMQ) so it runs       #
#  even while other FUSE files are in flux.                        #
#                                                                  #
#  Run with:                                                       #
#    cd FUSE && julia --project test/test_pedestal_nn_predictor.jl #
#= ============================================================ =#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using Statistics: mean, std

const __FUSE__ = abspath(joinpath(@__DIR__, ".."))
include(joinpath(__FUSE__, "src", "actors", "pedestal", "nn_predictor.jl"))

@testset "PedestalNN ensemble" begin
    onnx_dir = get(ENV, PEDESTAL_NN_ENV,
        normpath(joinpath(__FUSE__, "..", "pedestal-predictor-onnx", "onnx_models")))

    canonical_slug = "edensfit89"
    if !isfile(joinpath(onnx_dir, canonical_slug, "mse_encoder.onnx"))
        @warn "Skipping PedestalNN smoke test: ONNX bundle artefacts not found at $onnx_dir/$canonical_slug/"
        return
    end

    nn = load_pedestal_nn(; onnx_dir)

    @test length(nn.bundles) == 5
    @test sort(collect(keys(nn.bundles))) == sort(["edensfit89", "hmode_89", "te_ped_89", "ti_ped_89", "t_rot_ped_89"])
    @test length(nn.signal_names) == 32
    @test nn.signal_names[1] == "pohm"
    @test nn.signal_names[end] == "bt"
    @test length(nn.signal_means) == 32 == length(nn.signal_stds)

    # Per-bundle metadata sanity checks (values from manifest.json + target_norm.json)
    bedens = nn.bundles["edensfit89"]
    @test bedens.task === :regression
    @test bedens.target == "edens_ped"
    @test bedens.dataset_version == "v1_446"
    @test bedens.history_features == 446
    @test 2.0f0 < bedens.target_mean < 3.0f0     # ~2.58
    @test 1.0f0 < bedens.target_std  < 2.5f0     # ~1.61

    bte = nn.bundles["te_ped_89"]
    @test bte.task === :regression
    @test bte.dataset_version == "v2_458"
    @test bte.history_features == 458
    @test 0.3f0 < bte.target_mean < 0.7f0        # ~0.516 keV
    @test 0.2f0 < bte.target_std  < 0.6f0        # ~0.410 keV
    @test bte.target_units == "keV"

    bti = nn.bundles["ti_ped_89"]
    @test bti.task === :regression
    @test 0.6f0 < bti.target_mean < 1.2f0        # ~0.902 keV
    @test bti.target_units == "keV"

    bt_rot = nn.bundles["t_rot_ped_89"]
    @test bt_rot.task === :regression
    @test 10.0f0 < bt_rot.target_mean < 25.0f0   # ~17.19 krad/s
    @test bt_rot.target_units == "krad/s"

    bhmode = nn.bundles["hmode_89"]
    @test bhmode.task === :classification
    @test bhmode.dataset_version == "v1_446"
    @test bhmode.default_threshold == 0.5f0

    # Cached: second load returns the same object
    @test load_pedestal_nn(; onnx_dir) === nn

    # ── Run the ensemble with mean inputs (z-scored zeros everywhere) ──
    T = 200
    result = predict_pedestal(nn; T)

    @test length(result.te_ped) == T
    @test length(result.ti_ped) == T
    @test length(result.t_rot_ped) == T
    @test length(result.edens_ped) == T
    @test length(result.hmode_logit_seq) == T
    @test length(result.hmode_prob_seq) == T
    @test length(result.is_h_mode_seq) == T
    @test result.predictions_physical === result.edens_ped
    @test size(result.machine_embed) == (1, 512)
    @test size(result.aux_embed) == (1, 64)

    # Output statistics — for constant inputs the FPE produces a flat trace,
    # but it can still vary slightly via positional encodings; we just take
    # the per-trace mean and bound it to a physically plausible window.
    function trace_mean(x::AbstractVector)
        v = filter(!isnan, x)
        return isempty(v) ? NaN32 : Float32(sum(v) / length(v))
    end

    te_mean = trace_mean(result.te_ped)
    ti_mean = trace_mean(result.ti_ped)
    tr_mean = trace_mean(result.t_rot_ped)
    ne_mean = trace_mean(result.edens_ped)
    pH = result.hmode_prob

    @test 0.0f0 <= te_mean <= 5.0f0          # keV
    @test 0.0f0 <= ti_mean <= 5.0f0          # keV
    @test -200.0f0 <= tr_mean <= 200.0f0     # krad/s
    @test 0.5f0 <= ne_mean <= 10.0f0         # 1e19 m^-3
    @test 0.0f0 <= pH <= 1.0f0

    @info """
    predict_pedestal (all-zero / trained-mean inputs):
      ne_ped = $(round(ne_mean;  digits=3)) × 10^19 m^-3   [min=$(round(minimum(result.edens_ped);   digits=3)), max=$(round(maximum(result.edens_ped);   digits=3))]
      Te_ped = $(round(te_mean;  digits=3)) keV            [min=$(round(minimum(result.te_ped);     digits=3)), max=$(round(maximum(result.te_ped);     digits=3))]
      Ti_ped = $(round(ti_mean;  digits=3)) keV            [min=$(round(minimum(result.ti_ped);     digits=3)), max=$(round(maximum(result.ti_ped);     digits=3))]
      T_rot  = $(round(tr_mean;  digits=3)) krad/s         [min=$(round(minimum(result.t_rot_ped);  digits=3)), max=$(round(maximum(result.t_rot_ped);  digits=3))]
      P(H)   = $(round(pH;       digits=3))   logit=$(round(result.hmode_logit; digits=3))   is_h_mode=$(result.is_h_mode)
    """

    # Normalization helper round-trip: raw zeros z-scored should give -mean/std
    raw_zeros = zeros(Float32, 3, 32)
    znorm = normalized_signals(nn, raw_zeros)
    @test size(znorm) == (3, 32)
    @test all(isapprox.(znorm[1, :], -nn.signal_means ./ (nn.signal_stds .+ 1f-8); atol=1f-5))

    # Backward-compat alias (predict_density forwards to predict_pedestal)
    legacy = predict_density(nn; T=8)
    @test length(legacy.predictions_physical) == 8

    # ── ZMQ → FPE channel mapping table sanity ──
    # Every channel used by build_fpe_sequences_from_aux (in pedestal_actor.jl)
    # must exist in nn.signal_names with the expected index. Catches any future
    # bundle re-export that renames or reorders these channels.
    expected_indices = Dict(
        "pohm" => 1, "pinj" => 2,           # mapped from :zmq_Pohm (W) and :zmq_Pnbi (W → MW)
        "ip" => 30, "ipspr15v" => 31, "bt" => 32,
        "gasa_cal" => 25, "gasb_cal" => 26, "gasc_cal" => 27, "gasd_cal" => 28, "gase_cal" => 29,
        "f1a" => 5,  "f2a" => 6,  "f3a" => 7,  "f4a" => 8,  "f5a" => 9,
        "f6a" => 10, "f7a" => 11, "f8a" => 12, "f9a" => 13,
        "f1b" => 14, "f2b" => 15, "f3b" => 16, "f4b" => 17, "f5b" => 18,
        "f6b" => 19, "f7b" => 20, "f8b" => 21, "f9b" => 22,
        "ecoila" => 23, "ecoilb" => 24,    # mapped from :zmq_I_coil[1] (PCECOILA) and [4] (PCECOILB)
    )
    for (name, idx) in expected_indices
        @test fpe_signal_index(nn, name) == idx
    end
    # Channels still on training mean (pending GSLite-side :zmq_Tnbi / :zmq_Pech)
    for name in ("tinj", "ech_total")
        @test fpe_signal_index(nn, name) > 0
    end

    # Per-channel mean-fill convention: raw values equal to nn.signal_means
    # z-score to exactly zero (this is the fallback for unmapped channels in
    # build_fpe_sequences_from_aux).
    raw_means = repeat(reshape(nn.signal_means, 1, 32), 4, 1)
    @test all(abs.(normalized_signals(nn, raw_means)) .< 1f-5)

    # mean_normalized_history's aux_features must use the
    # (bzn=0, dis=-1, cov=0) sentinel from HistoryManager._init_from_raw,
    # not zeros (which would falsely claim disrupt-list coverage).
    _, _, aux_default = mean_normalized_history()
    @test aux_default == reshape(Float32[0f0, -1f0, 0f0], 1, 3)
    @test AUX_DEFAULT_SENTINEL == (0f0, -1f0, 0f0)

    # ── load_history_npz: round-trip via a synthetic NPZ ──
    # Builds an in-memory NPZ that matches HistoryManager.save_state's contract
    # (history_stats / history_masks / shot_ids / bzn/dis/cov scalars /
    # n_history_features), exercises both v2_458 and v1_446 widths, and asserts
    # shapes + z-score behavior + aux_features wiring. A real fixture from
    # `run_inference.py --save-history` lives behind FUSE_PEDESTAL_NN_HISTORY_NPZ;
    # if set, we additionally run predict_pedestal on it for a sanity sweep.
    import NPZ as _NPZ_test
    let tmpdir = mktempdir()
        for (n_features, label) in ((458, "v2_458"), (446, "v1_446"))
            n_slots_in = 50          # matches _MSE_HISTORY_LENGTH
            stats_raw = randn(Float32, n_slots_in, n_features)
            masks_in  = ones(Float32, n_slots_in)
            shots_in  = collect(Int64, 200_000:200_000+n_slots_in-1)
            bzn, dis, cov = 1234.5, -1.0, 0.0
            npz_path = joinpath(tmpdir, "synth_$(label).npz")
            _NPZ_test.npzwrite(npz_path, Dict(
                "history_stats"      => stats_raw,
                "history_masks"      => masks_in,
                "shot_ids"           => shots_in,
                "bzn_seconds"        => [bzn],
                "disrupt_seconds"    => [dis],
                "disrupt_coverage"   => [cov],
                "n_history_features" => Int64[n_features],
            ))

            # 1. Raw load (no norm params) — warns but should succeed and not z-score.
            hist, hmask, aux = @test_logs (:warn,) load_history_npz(npz_path)
            @test size(hist)  == (1, 50, 458)
            @test size(hmask) == (1, 50)
            @test size(aux)   == (1, 3)
            @test aux[1,1] ≈ Float32(bzn)
            @test aux[1,2] ≈ Float32(dis)
            @test aux[1,3] ≈ Float32(cov)
            # Width-padding: trailing channels [n_features+1:458] must be zero.
            if n_features < 458
                @test all(iszero, hist[1, :, n_features+1:end])
            end
            # Round-trip values for present channels match raw_stats (no z-score).
            @test all(isapprox.(hist[1, :, 1:n_features], stats_raw; atol=0))

            # 2. With matching norm-params JSON: z-scored stats land near 0 ± 1.
            μ = mean(stats_raw; dims=1) |> vec
            σ = max.(std(stats_raw; dims=1) |> vec, 1f-6)
            norm_json = joinpath(tmpdir, "norm_$(label).json")
            open(norm_json, "w") do io
                JSON.print(io, Dict("history_means"=>collect(μ), "history_stds"=>collect(σ)))
            end
            histz, _, _ = load_history_npz(npz_path; norm_params_path=norm_json)
            # z-scored values along the feature axis should have ~0 mean and ~1 std
            μ_z = vec(sum(histz[1, :, 1:n_features]; dims=1)) ./ size(histz, 2)
            @test maximum(abs, μ_z) < 1f-3
        end

        # 3. Slot truncation: 60 slots in → take latest 50, right-aligned.
        big = randn(Float32, 60, 458)
        big_masks = ones(Float32, 60)
        big_path = joinpath(tmpdir, "synth_60slots.npz")
        _NPZ_test.npzwrite(big_path, Dict(
            "history_stats"      => big,
            "history_masks"      => big_masks,
            "shot_ids"           => collect(Int64, 1:60),
            "bzn_seconds"        => [0.0],
            "disrupt_seconds"    => [-1.0],
            "disrupt_coverage"   => [0.0],
            "n_history_features" => Int64[458],
        ))
        hist60, _, _ = @test_logs (:warn,) load_history_npz(big_path)
        @test all(isapprox.(hist60[1, :, :], big[end-49:end, :]; atol=0))

        # 4. Slot zero-padding: 30 slots in → pad first 20 rows with zeros, mask=0.
        small = randn(Float32, 30, 458)
        small_path = joinpath(tmpdir, "synth_30slots.npz")
        _NPZ_test.npzwrite(small_path, Dict(
            "history_stats"      => small,
            "history_masks"      => ones(Float32, 30),
            "shot_ids"           => collect(Int64, 1:30),
            "bzn_seconds"        => [0.0],
            "disrupt_seconds"    => [-1.0],
            "disrupt_coverage"   => [0.0],
            "n_history_features" => Int64[458],
        ))
        hist30, hmask30, _ = @test_logs (:warn,) load_history_npz(small_path)
        @test all(iszero, hist30[1, 1:20, :])
        @test all(iszero, hmask30[1, 1:20])
        @test all(isone, hmask30[1, 21:50])
        @test all(isapprox.(hist30[1, 21:50, :], small; atol=0))
    end

    # ── Live history buffer round-trip (push_shot_history! → mse_history_from_aux) ──
    # Stand up a minimal `dd`-like that shares the `_aux::Dict{Symbol,Any}`
    # field accessed via `getfield(dd, :_aux)`. Avoids spinning up an IMAS.dd
    # for this skeleton test.
    mutable struct _MockDD
        _aux::Dict{Symbol,Any}
    end
    let mock = _MockDD(Dict{Symbol,Any}())
        # Empty buffer → mse_history_from_aux falls back to mean_normalized_history
        # (which warns about z-scoring? no — it doesn't z-score, it just returns zeros).
        h0, m0, a0 = mse_history_from_aux(mock, nn)
        @test size(h0) == (1, 50, 458)
        @test size(m0) == (1, 50)
        @test all(iszero, h0)
        @test all(isone, m0)             # mean fallback flags every slot as valid
        @test a0 == reshape(Float32[0f0, -1f0, 0f0], 1, 3)

        # Push 3 synthetic shots, each with distinct stats. With < 50 entries,
        # the buffer should right-align and zero-pad the earlier slots.
        s1 = randn(Float32, 458)
        s2 = randn(Float32, 458)
        s3 = randn(Float32, 458)
        push_shot_history!(mock, s1; shot_id=200_001, bzn_seconds=10.0)
        push_shot_history!(mock, s2; shot_id=200_002)
        push_shot_history!(mock, s3; shot_id=200_003,
                           disrupt_seconds=-1.0, disrupt_coverage=0.0)
        buf = mock._aux[:nn_history_buffer]
        @test buf.shot_ids == Int64[200_001, 200_002, 200_003]
        @test size(buf.history_stats) == (3, 458)
        @test buf.bzn_seconds == [10.0]              # set on first push, preserved through next two

        # mse_history_from_aux assembles the (1, 50, 458) tensor with the 3
        # latest rows right-aligned at slots 48-50 and zeros + mask=0 elsewhere.
        h, m, a = mse_history_from_aux(mock, nn)
        @test all(iszero, h[1, 1:47, :])
        @test all(iszero, m[1, 1:47])
        @test all(isone,  m[1, 48:50])
        @test all(isapprox.(h[1, 48, :], s1; atol=0))
        @test all(isapprox.(h[1, 49, :], s2; atol=0))
        @test all(isapprox.(h[1, 50, :], s3; atol=0))
        @test a[1,1] ≈ 10f0    # bzn from first push
        @test a[1,2] ≈ -1f0    # dis sentinel from third push
        @test a[1,3] ≈ 0f0     # cov sentinel from third push

        # Cap-at-50 behavior: push 60 more → buffer keeps only the latest 50.
        for k in 1:60
            push_shot_history!(mock, randn(Float32, 458); shot_id=300_000 + k)
        end
        buf2 = mock._aux[:nn_history_buffer]
        @test size(buf2.history_stats, 1) == 50
        @test buf2.shot_ids[end] == 300_060
        @test buf2.shot_ids[1]   == 300_011

        # Round-trip: save_history_npz → load_history_npz reproduces the buffer.
        let tmp = mktempdir()
            npz_out = joinpath(tmp, "rt.npz")
            save_history_npz(mock, npz_out)
            h_rt, m_rt, a_rt = load_history_npz(npz_out)
            h_buf, m_buf, a_buf = mse_history_from_aux(mock, nn)
            @test h_rt == h_buf
            @test m_rt == m_buf
            @test a_rt == a_buf
        end

        # compute_shot_stats stub returns zeros + warns; v1_446 vs v2_458 width.
        s_v2 = compute_shot_stats(mock; dataset_version="v2_458")
        s_v1 = compute_shot_stats(mock; dataset_version="v1_446")
        @test length(s_v2) == 458
        @test length(s_v1) == 446
        @test all(iszero, s_v2)
        @test_throws ErrorException compute_shot_stats(mock; dataset_version="bogus")
    end

    # Optional: real fixture from `run_inference.py --save-history`. Skipped if absent.
    # No norm_params_path is passed by default: none of the per-bundle
    # `normalization_params.json` files in onnx_models/ ship `history_means` /
    # `history_stds` keys, so the Python `HistoryManager` also forwards raw
    # stats to the encoder. Override with FUSE_PEDESTAL_NN_HISTORY_NORM_JSON if
    # a future bundle adds the per-feature norm vectors.
    history_npz = get(ENV, "FUSE_PEDESTAL_NN_HISTORY_NPZ", "")
    if !isempty(history_npz) && isfile(history_npz)
        @info "Running predict_pedestal on real history NPZ: $history_npz"
        norm_override = get(ENV, "FUSE_PEDESTAL_NN_HISTORY_NORM_JSON", "")
        norm_path = (!isempty(norm_override) && isfile(norm_override)) ? norm_override : nothing
        history = if norm_path === nothing
            @test_logs (:warn,) load_history_npz(history_npz)
        else
            load_history_npz(history_npz; norm_params_path=norm_path)
        end
        out = predict_pedestal(nn; history)
        for k in (:edens_ped, :te_ped, :ti_ped, :t_rot_ped)
            v = filter(isfinite, getfield(out, k))
            @test !isempty(v)
        end
    end

    # ── HuggingFace auto-fetch surface (offline-only checks) ──────────────
    @testset "download_pedestal_nn! manifest + cache-detection" begin
        files = pedestal_nn_files()
        # 1 top-level (manifest.json) + 4 regression bundles × 5 files + 1 classification × 4
        @test length(files) == 1 + 4 * 5 + 1 * 4
        @test "manifest.json" in files
        @test "edensfit89/mse_encoder.onnx" in files
        @test "edensfit89/fpe_encoder.onnx" in files
        @test "edensfit89/target_norm.json" in files
        # hmode_89 is classification — no target_norm.json
        @test "hmode_89/mse_encoder.onnx" in files
        @test !("hmode_89/target_norm.json" in files)
        # All 5 canonical slugs present
        for slug in _CANONICAL_SLUGS
            @test any(startswith(f, slug * "/") for f in files)
        end

        # pedestal_nn_dir_complete: false on empty dir, true on a fake-populated one
        mktempdir() do empty_tmp
            @test pedestal_nn_dir_complete(empty_tmp) === false
        end
        mktempdir() do fake_tmp
            for rel in files
                p = joinpath(fake_tmp, rel)
                mkpath(dirname(p))
                write(p, "x")  # non-empty stub
            end
            @test pedestal_nn_dir_complete(fake_tmp) === true
            # zero-byte file should make it incomplete again
            zero_target = joinpath(fake_tmp, files[end])
            write(zero_target, "")
            @test pedestal_nn_dir_complete(fake_tmp) === false
        end

        # Auto-fetch hook: explicit opt-out (`0`/`false`/`no`/`off`) skips the
        # download even on an empty dir. No network call expected.
        for v in ("0", "false", "no", "off", "FALSE", "Off")
            mktempdir() do tmp
                old = get(ENV, PEDESTAL_NN_AUTOFETCH_ENV, nothing)
                ENV[PEDESTAL_NN_AUTOFETCH_ENV] = v
                try
                    @test _maybe_autofetch!(tmp) === false
                finally
                    old === nothing ? delete!(ENV, PEDESTAL_NN_AUTOFETCH_ENV) :
                                      (ENV[PEDESTAL_NN_AUTOFETCH_ENV] = old)
                end
            end
        end

        # Auto-fetch hook: when the dir is *already complete*, no fetch is
        # attempted regardless of the env var (idempotency). This is the
        # branch that protects against repeated downloads on every container
        # restart and against `using FUSE` re-fetching when files exist.
        for v in ("", "1", "true", "0", "false")
            mktempdir() do tmp
                # Stub a complete bundle dir.
                for rel in pedestal_nn_files()
                    p = joinpath(tmp, rel)
                    mkpath(dirname(p))
                    write(p, "x")
                end
                @assert pedestal_nn_dir_complete(tmp)
                old = get(ENV, PEDESTAL_NN_AUTOFETCH_ENV, nothing)
                if isempty(v)
                    delete!(ENV, PEDESTAL_NN_AUTOFETCH_ENV)
                else
                    ENV[PEDESTAL_NN_AUTOFETCH_ENV] = v
                end
                try
                    @test _maybe_autofetch!(tmp) === false
                finally
                    old === nothing ? delete!(ENV, PEDESTAL_NN_AUTOFETCH_ENV) :
                                      (ENV[PEDESTAL_NN_AUTOFETCH_ENV] = old)
                end
            end
        end

        # Optional: real network round-trip, gated by env var (CI doesn't run).
        # Set FUSE_PEDESTAL_NN_DOWNLOAD_TEST=1 to fetch a single small JSON
        # from HuggingFace and verify the URL/header pipeline works.
        if get(ENV, "FUSE_PEDESTAL_NN_DOWNLOAD_TEST", "0") == "1"
            mktempdir() do tmp
                # download just the cheapest file (manifest.json, ~few KB) by
                # temporarily monkey-patching the file manifest at the call site
                # via a stripped-down call.
                base = "https://huggingface.co/$(PEDESTAL_NN_HF_REPO_DEFAULT)/resolve/$(PEDESTAL_NN_HF_REVISION_DEFAULT)"
                dest = joinpath(tmp, "manifest.json")
                Downloads.download("$base/manifest.json", dest)
                @test isfile(dest)
                @test filesize(dest) > 0
                m = JSON.parsefile(dest)
                @test haskey(m, "bundles")
                @test "edensfit89" in keys(m["bundles"])
            end
        end
    end
end
