#= =========================================================== =#
#  Smoke test for the PedestalDensityNN wrapper.                  #
#                                                                 #
#  Activates the FUSE project and exercises the nn_predictor.jl   #
#  module in isolation (no actors, no dd, no ZMQ) so it runs      #
#  even while other FUSE files are in flux.                       #
#                                                                 #
#  Run with:                                                      #
#    cd FUSE && julia --project test/test_pedestal_nn_predictor.jl #
#= =========================================================== =#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test

const __FUSE__ = abspath(joinpath(@__DIR__, ".."))
include(joinpath(__FUSE__, "src", "actors", "pedestal", "nn_predictor.jl"))

@testset "PedestalDensityNN" begin
    onnx_dir = get(ENV, PEDESTAL_NN_ENV,
        normpath(joinpath(__FUSE__, "..", "pedestal-predictor-onnx", "onnx_models")))

    if !isfile(joinpath(onnx_dir, "mse_encoder.onnx"))
        @warn "Skipping PedestalDensityNN smoke test: ONNX artefacts not found at $onnx_dir"
        return
    end

    nn = load_pedestal_nn(; onnx_dir)
    @test length(nn.signal_names) == 32
    @test nn.signal_names[1] == "pohm"
    @test nn.signal_names[end] == "bt"
    @test length(nn.signal_means) == 32
    @test length(nn.signal_stds) == 32
    @test 2.0f0 < nn.target_mean < 3.0f0    # edensfit89 training mean ~ 2.58
    @test 1.0f0 < nn.target_std < 2.5f0     #                      std ~ 1.61

    # Cached: second load returns the same object
    @test load_pedestal_nn(; onnx_dir) === nn

    result = predict_density(nn; T=200)

    @test size(result.machine_embed) == (1, 512)
    @test size(result.aux_embed) == (1, 64)
    @test length(result.predictions_physical) == 200
    @test length(result.predictions_norm) == 200

    m = sum(result.predictions_physical) / length(result.predictions_physical)
    @test 0.5f0 <= m <= 10.0f0   # physical plausibility window from the Python pipeline
    @info "predict_density (all-zero / trained-mean inputs): mean=$(round(m; digits=3))  min=$(round(minimum(result.predictions_physical); digits=3))  max=$(round(maximum(result.predictions_physical); digits=3))  x 10^19 m^-3"

    # Normalization helper round-trip: raw zeros z-scored should give -mean/std
    raw = zeros(Float32, 3, 32)
    znorm = normalized_signals(nn, raw)
    @test size(znorm) == (3, 32)
    @test all(isapprox.(znorm[1, :], -nn.signal_means ./ (nn.signal_stds .+ 1f-8); atol=1f-5))
end
