#= =========================================================== =#
#  PedestalNNPredictor                                            #
#                                                                 #
#  Edensfit89 two-stage ONNX pipeline for pedestal electron       #
#  density prediction on DIII-D (Harris et al., FUSIONDL AESP).   #
#                                                                 #
#    MSE encoder : (history_stats (1,50,446) +                    #
#                   history_masks (1,50) +                        #
#                   aux_features  (1,3)) -> machine_embed (1,512) #
#                                        + aux_embed     (1, 64)  #
#    FPE encoder : (sequences    (1, T, 32) +                     #
#                   machine_embed,                                #
#                   aux_embed,                                    #
#                   signal_masks (1, 32),                         #
#                   padding_mask (1, T)) -> predictions (1, T)    #
#                                                                 #
#  Physical density (10^19 m^-3) = pred * target_std + target_mean #
#                                                                 #
#  The ONNX files and normalization JSONs live *outside* FUSE      #
#  (~148 MB of weights). Their directory is resolved from, in      #
#  order:                                                          #
#    (1) the `onnx_dir` kwarg passed to `load_pedestal_nn`,        #
#    (2) the environment variable `FUSE_PEDESTAL_NN_DIR`,          #
#    (3) `<pedestal-predictor-onnx repo>/onnx_models/` alongside   #
#        the FUSE development tree.                                #
#                                                                 #
#  Runtime: ONNXRunTime.jl (Microsoft onnxruntime C library).     #
#  ONNXNaiveNASflux cannot load these PyTorch opset-17 graphs     #
#  because they use Cast, ConstantOfShape, ScatterND, masked      #
#  Unsqueeze patterns, etc.                                        #
#= =========================================================== =#

import ONNXRunTime
import JSON

const PEDESTAL_NN_ENV = "FUSE_PEDESTAL_NN_DIR"

"""
    PedestalDensityNN

Loaded edensfit89 pipeline: the two ONNX sessions plus the normalization
parameters needed to map physical inputs <-> z-scored network I/O.

Fields:
- `mse`               MSE-encoder inference session
- `fpe`               FPE-encoder inference session
- `signal_names`      ordered names of the 32 FPE channels
- `signal_means`      per-channel FPE training means  (Vector{Float32}, 32)
- `signal_stds`       per-channel FPE training stdevs (Vector{Float32}, 32)
- `target_mean`       pedestal density training mean (10^19 m^-3)
- `target_std`        pedestal density training std  (10^19 m^-3)
- `onnx_dir`          directory the artefacts were loaded from (for reporting)
"""
struct PedestalDensityNN
    mse::ONNXRunTime.InferenceSession
    fpe::ONNXRunTime.InferenceSession
    signal_names::Vector{String}
    signal_means::Vector{Float32}
    signal_stds::Vector{Float32}
    target_mean::Float32
    target_std::Float32
    onnx_dir::String
end

"""
    resolve_pedestal_nn_dir(; onnx_dir=nothing) -> String

Resolve the directory that holds `mse_encoder.onnx`, `fpe_encoder.onnx`,
`normalization_params.json`, and `dataset_metadata.json`.
"""
function resolve_pedestal_nn_dir(; onnx_dir::Union{Nothing,AbstractString}=nothing)
    if onnx_dir !== nothing
        return String(onnx_dir)
    end
    if haskey(ENV, PEDESTAL_NN_ENV)
        return ENV[PEDESTAL_NN_ENV]
    end
    # Default: sibling checkout of pedestal-predictor-onnx next to FUSE
    default = normpath(joinpath(__FUSE__, "..", "pedestal-predictor-onnx", "onnx_models"))
    return default
end

"""
    load_pedestal_nn(; onnx_dir=nothing) -> PedestalDensityNN

Load the edensfit89 ONNX pipeline and the companion JSON metadata. Cached
per-process: subsequent calls with the same directory return the same object.

Set `ENV["FUSE_PEDESTAL_NN_DIR"]` or pass `onnx_dir` to point at a non-default
location. The directory must contain:

    mse_encoder.onnx
    fpe_encoder.onnx
    normalization_params.json
    dataset_metadata.json
"""
function load_pedestal_nn(; onnx_dir::Union{Nothing,AbstractString}=nothing)
    dir = resolve_pedestal_nn_dir(; onnx_dir)
    return _cached_pedestal_nn(dir)
end

const _PEDESTAL_NN_CACHE = Dict{String,PedestalDensityNN}()

function _cached_pedestal_nn(dir::AbstractString)
    key = abspath(dir)
    m = get(_PEDESTAL_NN_CACHE, key, nothing)
    m === nothing || return m
    isdir(key) || error("PedestalDensityNN: onnx directory does not exist: $key\n" *
                        "Set ENV[\"$PEDESTAL_NN_ENV\"] or pass onnx_dir=... to load_pedestal_nn.")

    mse_path = joinpath(key, "mse_encoder.onnx")
    fpe_path = joinpath(key, "fpe_encoder.onnx")
    norm_path = joinpath(key, "normalization_params.json")
    meta_path = joinpath(key, "dataset_metadata.json")
    for p in (mse_path, fpe_path, norm_path, meta_path)
        isfile(p) || error("PedestalDensityNN: missing artefact: $p")
    end

    mse = ONNXRunTime.load_inference(mse_path)
    fpe = ONNXRunTime.load_inference(fpe_path)

    norm = JSON.parsefile(norm_path)
    signal_names = String.(norm["signal_names"])
    signal_means = Float32.(norm["means"])
    signal_stds = Float32.(norm["stds"])
    length(signal_means) == 32 || error("expected 32 FPE means, got $(length(signal_means))")
    length(signal_stds) == 32 || error("expected 32 FPE stds, got $(length(signal_stds))")

    meta = JSON.parsefile(meta_path)
    target_mean = Float32(meta["target_mean"])
    target_std = Float32(meta["target_std"])

    m = PedestalDensityNN(mse, fpe, signal_names, signal_means, signal_stds, target_mean, target_std, key)
    _PEDESTAL_NN_CACHE[key] = m
    return m
end

"""
    fpe_signal_index(nn::PedestalDensityNN, name::AbstractString) -> Int

Return the FPE channel index for `name` (e.g. "pinj", "bt"), or 0 if absent.
"""
function fpe_signal_index(nn::PedestalDensityNN, name::AbstractString)
    idx = findfirst(==(String(name)), nn.signal_names)
    idx === nothing ? 0 : idx
end

"""
    normalized_signals(nn::PedestalDensityNN, raw::AbstractMatrix) -> Matrix{Float32}

Z-score a `(T, 32)` matrix of raw-physical FPE signals using the training
per-channel mean/std stored in `nn`.
"""
function normalized_signals(nn::PedestalDensityNN, raw::AbstractMatrix)
    size(raw, 2) == length(nn.signal_means) ||
        error("normalized_signals: expected $(length(nn.signal_means)) channels, got $(size(raw,2))")
    T = size(raw, 1)
    out = Matrix{Float32}(undef, T, length(nn.signal_means))
    @inbounds for c in axes(raw, 2), t in axes(raw, 1)
        out[t, c] = (Float32(raw[t, c]) - nn.signal_means[c]) / (nn.signal_stds[c] + 1f-8)
    end
    return out
end

"""
    mean_normalized_history() -> (hist, hmask, aux)

Produce the "neutral" MSE inputs corresponding to the normalization means
(z-score zeros). Use this when no machine-state history file is available.
Shapes match the ONNX graph: `(1, 50, 446)`, `(1, 50)`, `(1, 3)`.
"""
function mean_normalized_history()
    hist  = zeros(Float32, 1, 50, 446)
    hmask = ones(Float32, 1, 50)            # all 50 history slots valid
    aux   = zeros(Float32, 1, 3)            # 0 = z-scored mean for bzn/disrupt aux
    return hist, hmask, aux
end

"""
    load_history_npz(path) -> (hist, hmask, aux)

Load an MSE history tensor previously saved by the Python pipeline's
`--save-history` flag. Returns arrays shaped to match the ONNX graph inputs.
Not yet wired — stubbed for the next milestone.
"""
function load_history_npz(path::AbstractString)
    error("load_history_npz: not yet implemented. Pass mean_normalized_history() until MSE history loader is ready.")
end

"""
    predict_density(nn::PedestalDensityNN;
                    sequences::Union{Nothing,AbstractMatrix}=nothing,
                    T::Integer=200,
                    prenormalized::Bool=false,
                    signal_mask::Union{Nothing,AbstractVector}=nothing,
                    padding_mask::Union{Nothing,AbstractVector}=nothing,
                    history::Union{Nothing,Tuple}=nothing)
        -> NamedTuple

Run the two-stage MSE + FPE pipeline and return the predicted pedestal
electron density trace.

Arguments:
- `sequences`  : `(T, 32)` matrix of FPE actuator signals. If `nothing`,
                 uses an all-zeros (= z-scored mean) sequence of length `T`.
- `prenormalized` : set `true` when `sequences` is already z-scored. Default
                 `false` means raw physical values are z-scored inside.
- `signal_mask` : length-32 Float32 vector, 1.0 where the channel is present.
                 Defaults to all-ones.
- `padding_mask` : length-T Float32 vector, 1.0 for valid timesteps. Defaults
                  to all-ones.
- `history` : `(hist, hmask, aux)` triple of MSE inputs. Defaults to
              `mean_normalized_history()`.

Returns NamedTuple with:
- `predictions_physical::Vector{Float32}` — density in 10^19 m^-3, length T
- `predictions_norm::Vector{Float32}` — z-scored model output, length T
- `machine_embed::Matrix{Float32}` — (1, 512)
- `aux_embed::Matrix{Float32}` — (1, 64)
"""
function predict_density(nn::PedestalDensityNN;
    sequences::Union{Nothing,AbstractMatrix}=nothing,
    T::Integer=200,
    prenormalized::Bool=false,
    signal_mask::Union{Nothing,AbstractVector}=nothing,
    padding_mask::Union{Nothing,AbstractVector}=nothing,
    history::Union{Nothing,Tuple}=nothing)

    # ── MSE inputs ────────────────────────────────────────────────
    if history === nothing
        hist, hmask, aux = mean_normalized_history()
    else
        hist, hmask, aux = history
    end
    size(hist) == (1, 50, 446) || error("history_stats must be (1, 50, 446), got $(size(hist))")
    size(hmask) == (1, 50) || error("history_masks must be (1, 50), got $(size(hmask))")
    size(aux) == (1, 3) || error("aux_features must be (1, 3), got $(size(aux))")

    mse_out = nn.mse(Dict(
        "history_stats" => Float32.(hist),
        "history_masks" => Float32.(hmask),
        "aux_features" => Float32.(aux)))
    machine_embed = mse_out["machine_embed"]::AbstractMatrix
    aux_embed = mse_out["aux_embed"]::AbstractMatrix

    # ── FPE inputs ────────────────────────────────────────────────
    if sequences === nothing
        # z-scored zero = per-channel training mean
        sequences_norm = zeros(Float32, T, 32)
    else
        size(sequences, 2) == 32 || error("sequences must have 32 channels, got $(size(sequences,2))")
        T = size(sequences, 1)
        sequences_norm = prenormalized ? Matrix{Float32}(sequences) :
                         normalized_signals(nn, sequences)
    end
    smask = signal_mask === nothing ? ones(Float32, 32) : Float32.(signal_mask)
    length(smask) == 32 || error("signal_mask must have length 32")
    pmask = padding_mask === nothing ? ones(Float32, T) : Float32.(padding_mask)
    length(pmask) == T || error("padding_mask length ($(length(pmask))) must equal T ($T)")

    seq3 = reshape(sequences_norm, 1, T, 32)
    smask2 = reshape(smask, 1, 32)
    pmask2 = reshape(pmask, 1, T)

    fpe_out = nn.fpe(Dict(
        "sequences" => seq3,
        "machine_embed" => machine_embed,
        "aux_embed" => aux_embed,
        "signal_masks" => smask2,
        "padding_mask" => pmask2))
    pred_norm = vec(fpe_out["predictions"]::AbstractMatrix)
    pred_phys = @. pred_norm * nn.target_std + nn.target_mean

    return (;
        predictions_physical=Vector{Float32}(pred_phys),
        predictions_norm=Vector{Float32}(pred_norm),
        machine_embed=Matrix{Float32}(machine_embed),
        aux_embed=Matrix{Float32}(aux_embed))
end
