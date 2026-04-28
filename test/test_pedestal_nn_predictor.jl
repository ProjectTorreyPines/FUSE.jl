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
end
