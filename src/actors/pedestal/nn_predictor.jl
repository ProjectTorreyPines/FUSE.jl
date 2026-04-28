#= ============================================================================ =#
#  PedestalNN — DIII-D pedestal NN ensemble (PedestalPredictor v2)              #
#                                                                                #
#  Five-bundle ONNX ensemble published at                                        #
#    https://huggingface.co/SCS-Lab/pedestal-predictor-onnx                      #
#                                                                                #
#    edensfit89   regression  -> edens_ped   (10^19 m^-3, v1_446)                #
#    hmode_89     classification -> hmode    (probability,  v1_446)              #
#    te_ped_89    regression  -> te_ped      (keV,          v2_458)              #
#    ti_ped_89    regression  -> ti_ped      (keV,          v2_458)              #
#    t_rot_ped_89 regression  -> t_rot_ped   (krad/s,       v2_458)              #
#                                                                                #
#  Each bundle has its own MSE + FPE encoder pair sharing the same architecture, #
#  the same FPE input contract (32 actuators), and (per docs / verified at load) #
#  the same per-channel FPE normalization. The MSE history width differs:        #
#  v1 bundles want 446-wide history, v2 bundles want 458-wide. We keep the       #
#  superset (458) internally and slice to 446 for v1 bundles.                    #
#                                                                                #
#  Runtime: ONNXRunTime.jl (Microsoft onnxruntime C library). The PyTorch        #
#  opset-17 graphs use Cast / ConstantOfShape / ScatterND / masked-Unsqueeze     #
#  patterns that ONNXNaiveNASflux cannot lower.                                  #
#                                                                                #
#  Artefacts live *outside* FUSE (~700 MB of weights). Their root directory is   #
#  resolved from, in order:                                                      #
#    (1) the `onnx_dir` kwarg passed to `load_pedestal_nn`,                      #
#    (2) the environment variable `FUSE_PEDESTAL_NN_DIR`,                        #
#    (3) `<pedestal-predictor-onnx repo>/onnx_models/` alongside FUSE.           #
#= ============================================================================ =#

import ONNXRunTime
import JSON
import OrderedCollections: OrderedDict

const PEDESTAL_NN_ENV = "FUSE_PEDESTAL_NN_DIR"

const _CANONICAL_SLUGS = ("edensfit89", "hmode_89", "te_ped_89", "ti_ped_89", "t_rot_ped_89")
const _DATASET_WIDTH = Dict("v1_446" => 446, "v2_458" => 458)
const _MSE_HISTORY_LENGTH = 50
const _FPE_SIGNAL_DIM = 32
const _AUX_DIM = 3
const _MAX_HISTORY_WIDTH = 458    # v2_458 superset; v1_446 bundles get [..., 1:446]

"""
    PedestalNNBundle

One ONNX bundle (MSE + FPE encoder pair) plus its sidecar metadata.
"""
struct PedestalNNBundle
    slug::String
    task::Symbol                       # :regression | :classification
    target::String                     # te_ped, ti_ped, t_rot_ped, edens_ped, hmode
    dataset_version::String            # v1_446 | v2_458
    history_features::Int
    signal_names::Vector{String}
    means::Vector{Float32}
    stds::Vector{Float32}
    target_mean::Float32               # NaN32 for classification
    target_std::Float32                # NaN32 for classification
    target_units::String               # "" for classification
    default_threshold::Float32         # 0.5 for classification, NaN32 for regression
    mse::ONNXRunTime.InferenceSession
    fpe::ONNXRunTime.InferenceSession
end

is_regression(b::PedestalNNBundle) = b.task === :regression
is_classification(b::PedestalNNBundle) = b.task === :classification

"""
    PedestalNN

Loaded PedestalPredictor ONNX ensemble. Exposes one inference entry point —
[`predict_pedestal`](@ref) — that runs all loaded bundles on a shared set of
inputs and returns physical-unit predictions for ne_ped, Te_ped, Ti_ped,
T_rot_ped, and the H-mode probability.

Fields:
- `bundles`      : ordered slug -> [`PedestalNNBundle`](@ref)
- `signal_names` : 32 FPE channel names shared by every bundle
- `signal_means` : per-channel FPE training means  (Float32 × 32)
- `signal_stds`  : per-channel FPE training stdevs (Float32 × 32)
- `onnx_dir`     : directory the artefacts were loaded from (for reporting)
"""
struct PedestalNN
    bundles::OrderedDict{String,PedestalNNBundle}
    signal_names::Vector{String}
    signal_means::Vector{Float32}
    signal_stds::Vector{Float32}
    onnx_dir::String
end

# Backward-compat alias for any consumer that still references the old type
const PedestalDensityNN = PedestalNN

# ── Path resolution & loading ───────────────────────────────────────────────

"""
    resolve_pedestal_nn_dir(; onnx_dir=nothing) -> String

Resolve the directory that holds the per-bundle ONNX subdirs
(`edensfit89/`, `hmode_89/`, `te_ped_89/`, `ti_ped_89/`, `t_rot_ped_89/`).
"""
function resolve_pedestal_nn_dir(; onnx_dir::Union{Nothing,AbstractString}=nothing)
    if onnx_dir !== nothing
        return String(onnx_dir)
    end
    if haskey(ENV, PEDESTAL_NN_ENV)
        return ENV[PEDESTAL_NN_ENV]
    end
    return normpath(joinpath(__FUSE__, "..", "pedestal-predictor-onnx", "onnx_models"))
end

"""
    load_pedestal_nn(; onnx_dir=nothing, only=nothing) -> PedestalNN

Load every available PedestalPredictor ONNX bundle from `onnx_dir` (resolved
via [`resolve_pedestal_nn_dir`](@ref)). Cached per-process: subsequent calls
with the same directory return the same object.

Pass `only=["edensfit89", "hmode_89"]` to restrict to a subset of bundles
(the others' fields in [`predict_pedestal`](@ref)'s NamedTuple will be NaN).

Set `ENV["FUSE_PEDESTAL_NN_DIR"]` or pass `onnx_dir=...` to point at a
non-default location. The directory must contain one subdirectory per
loaded slug, each holding `mse_encoder.onnx`, `fpe_encoder.onnx`,
`model_config.json`, `normalization_params.json`, and (for regression
bundles) `target_norm.json`.
"""
function load_pedestal_nn(; onnx_dir::Union{Nothing,AbstractString}=nothing,
                            only::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    dir = resolve_pedestal_nn_dir(; onnx_dir)
    return _cached_pedestal_nn(dir, only)
end

const _PEDESTAL_NN_CACHE = Dict{Tuple{String,Vector{String}},PedestalNN}()

function _cached_pedestal_nn(dir::AbstractString, only)
    keep = only === nothing ? collect(_CANONICAL_SLUGS) : sort(String.(only))
    cache_key = (abspath(dir), keep)
    haskey(_PEDESTAL_NN_CACHE, cache_key) && return _PEDESTAL_NN_CACHE[cache_key]
    isdir(cache_key[1]) || error("PedestalNN: onnx directory does not exist: $(cache_key[1])\n" *
                                 "Set ENV[\"$PEDESTAL_NN_ENV\"] or pass onnx_dir=... to load_pedestal_nn.")

    bundles = OrderedDict{String,PedestalNNBundle}()
    for slug in _CANONICAL_SLUGS
        only === nothing || (slug in only) || continue
        bdir = joinpath(cache_key[1], slug)
        isdir(bdir) || (@warn "PedestalNN: bundle directory missing, skipping: $bdir"; continue)
        bundles[slug] = _load_bundle(slug, bdir)
    end
    isempty(bundles) && error("PedestalNN: no bundles loaded from $(cache_key[1]). " *
                              "Expected one or more of: $(_CANONICAL_SLUGS)")

    ref = first(values(bundles))
    for b in values(bundles)
        b.signal_names == ref.signal_names ||
            error("PedestalNN: $(b.slug) FPE signal_names differ from $(ref.slug) — cannot run a joint ensemble.")
        maximum(abs.(b.means .- ref.means)) < 1f-5 ||
            @warn "PedestalNN: $(b.slug) FPE means differ from $(ref.slug) by up to $(maximum(abs.(b.means .- ref.means)))"
        maximum(abs.(b.stds .- ref.stds)) < 1f-5 ||
            @warn "PedestalNN: $(b.slug) FPE stds differ from $(ref.slug) by up to $(maximum(abs.(b.stds .- ref.stds)))"
    end

    nn = PedestalNN(bundles, ref.signal_names, ref.means, ref.stds, cache_key[1])
    _PEDESTAL_NN_CACHE[cache_key] = nn
    return nn
end

function _load_bundle(slug::AbstractString, bdir::AbstractString)
    cfg_path = joinpath(bdir, "model_config.json")
    norm_path = joinpath(bdir, "normalization_params.json")
    mse_path = joinpath(bdir, "mse_encoder.onnx")
    fpe_path = joinpath(bdir, "fpe_encoder.onnx")
    for p in (cfg_path, norm_path, mse_path, fpe_path)
        isfile(p) || error("PedestalNN $slug: missing artefact: $p")
    end

    cfg = JSON.parsefile(cfg_path)
    norm = JSON.parsefile(norm_path)

    task_str = String(cfg["task_type"])
    task = task_str == "regression" ? :regression :
           task_str == "classification" ? :classification :
           error("$slug: invalid task_type=$(task_str)")

    dsv = String(cfg["dataset_version"])
    haskey(_DATASET_WIDTH, dsv) || error("$slug: invalid dataset_version=$(dsv)")
    expected_width = _DATASET_WIDTH[dsv]
    Int(cfg["history_features"]) == expected_width ||
        error("$slug: model_config.json::history_features $(cfg["history_features"]) " *
              "!= expected $expected_width for dataset_version=$dsv")

    signal_names = String.(norm["signal_names"])
    means = Float32.(norm["means"])
    stds = Float32.(norm["stds"])
    length(means) == _FPE_SIGNAL_DIM ||
        error("$slug: expected $_FPE_SIGNAL_DIM FPE means, got $(length(means))")
    length(stds) == _FPE_SIGNAL_DIM ||
        error("$slug: expected $_FPE_SIGNAL_DIM FPE stds, got $(length(stds))")

    target_mean = NaN32
    target_std = NaN32
    target_units = ""
    default_threshold = NaN32
    if task === :regression
        tmean = get(cfg, "target_mean", nothing)
        tstd = get(cfg, "target_std", nothing)
        if tmean === nothing || tstd === nothing
            tn_path = joinpath(bdir, "target_norm.json")
            isfile(tn_path) || error("$slug: regression bundle missing target_mean/target_std and no target_norm.json")
            tn = JSON.parsefile(tn_path)
            tmean = something(tmean, tn["target_mean"])
            tstd = something(tstd, tn["target_std"])
        end
        target_mean = Float32(tmean)
        target_std = Float32(tstd)
        target_units = String(get(cfg, "target_units", ""))
    else
        default_threshold = Float32(get(cfg, "default_threshold", 0.5))
    end

    mse = ONNXRunTime.load_inference(mse_path)
    fpe = ONNXRunTime.load_inference(fpe_path)

    return PedestalNNBundle(
        String(slug), task, String(cfg["target_name"]), dsv, expected_width,
        signal_names, means, stds, target_mean, target_std, target_units, default_threshold,
        mse, fpe)
end

# ── Normalization helpers ───────────────────────────────────────────────────

"""
    fpe_signal_index(nn::PedestalNN, name::AbstractString) -> Int

Return the FPE channel index for `name` (e.g. "pinj", "bt"), or 0 if absent.
"""
function fpe_signal_index(nn::PedestalNN, name::AbstractString)
    idx = findfirst(==(String(name)), nn.signal_names)
    return idx === nothing ? 0 : idx
end

"""
    normalized_signals(nn::PedestalNN, raw::AbstractMatrix) -> Matrix{Float32}

Z-score a `(T, 32)` matrix of raw-physical FPE signals using the training
per-channel mean/std stored in `nn`.
"""
function normalized_signals(nn::PedestalNN, raw::AbstractMatrix)
    size(raw, 2) == _FPE_SIGNAL_DIM ||
        error("normalized_signals: expected $_FPE_SIGNAL_DIM channels, got $(size(raw,2))")
    T = size(raw, 1)
    out = Matrix{Float32}(undef, T, _FPE_SIGNAL_DIM)
    @inbounds for c in 1:_FPE_SIGNAL_DIM, t in 1:T
        out[t, c] = (Float32(raw[t, c]) - nn.signal_means[c]) / (nn.signal_stds[c] + 1f-8)
    end
    return out
end

"""
    mean_normalized_history() -> (hist, hmask, aux)

Produce the "neutral" MSE inputs corresponding to the normalization means
(z-score zeros). Use this when no machine-state history file is available.
Shapes match the v2 ONNX graph: `(1, 50, 458)`, `(1, 50)`, `(1, 3)`.
v1 bundles internally see `[..., 1:446]`.
"""
function mean_normalized_history()
    hist = zeros(Float32, 1, _MSE_HISTORY_LENGTH, _MAX_HISTORY_WIDTH)
    hmask = ones(Float32, 1, _MSE_HISTORY_LENGTH)
    aux = zeros(Float32, 1, _AUX_DIM)
    return hist, hmask, aux
end

"""
    load_history_npz(path) -> (hist, hmask, aux)

Load an MSE history tensor previously saved by the Python pipeline's
`--save-history` flag. Returns arrays shaped to match the ONNX graph inputs.
Not yet wired — stubbed for the next milestone.
"""
function load_history_npz(::AbstractString)
    error("load_history_npz: not yet implemented. Pass mean_normalized_history() until the MSE history loader is ready.")
end

# ── Inference ───────────────────────────────────────────────────────────────

"""
    predict_pedestal(nn::PedestalNN; sequences=nothing, T=200, prenormalized=false,
                     signal_mask=nothing, padding_mask=nothing,
                     history=nothing, h_mode_threshold=nothing) -> NamedTuple

Run every loaded bundle on a shared set of inputs and return a NamedTuple of
*physical-unit* predictions plus advanced raw outputs. Defaults to "z-scored
zero" inputs (= per-channel training means) for both the MSE history and the
FPE actuator sequences.

Arguments:
- `sequences`     : `(T, 32)` matrix of FPE actuator signals. If `nothing`,
                    uses an all-zeros (= z-scored mean) sequence of length `T`.
- `prenormalized` : set `true` when `sequences` is already z-scored. Default
                    `false` means raw physical values are z-scored inside.
- `signal_mask`   : length-32 vector, 1.0 where the channel is present.
                    Defaults to all-ones.
- `padding_mask`  : length-T vector, 1.0 for valid timesteps. Defaults to
                    all-ones.
- `history`       : `(hist, hmask, aux)` triple of MSE inputs (full v2_458
                    width). Defaults to [`mean_normalized_history`](@ref).
- `h_mode_threshold` : override for the H-mode classifier threshold (default
                    `0.5` from `hmode_89/model_config.json::default_threshold`).

Returns a NamedTuple with:
- `te_ped`, `ti_ped`, `t_rot_ped`, `edens_ped` :: Vector{Float32} (length T)
   — denormalized physical-unit traces (filled with `NaN32` if the bundle
   wasn't loaded; units in the `units` dict).
- `hmode_logit_seq`, `hmode_prob_seq` :: Vector{Float32} (length T)
   — per-timestep raw logit and sigmoid probability.
- `is_h_mode_seq` :: Vector{Bool} (length T) — `hmode_prob_seq .>= threshold`.
- `hmode_logit`, `hmode_prob` :: Float32 — scalar (mean over T) summary used
   by `ActorPedestal` for L/H mode dispatch.
- `is_h_mode` :: Bool — `hmode_prob >= threshold`.
- `predictions_physical`, `predictions_norm` :: Vector{Float32}
   — backward-compat aliases for the density-only predictor (= `edens_ped`
   and the corresponding pre-denormalization output).
- `raw` :: Dict{String,Vector{Float32}} — per-slug pre-denorm / pre-sigmoid
   ONNX outputs, keyed by bundle slug.
- `units` :: Dict{Symbol,String} — physical-unit labels by output name.
- `machine_embed`, `aux_embed` :: Matrix{Float32} — kept for backward
   compatibility (taken from any bundle's MSE since they're independently
   trained but the same shape `(1, 512)` / `(1, 64)`).
"""
function predict_pedestal(nn::PedestalNN;
        sequences::Union{Nothing,AbstractMatrix}=nothing,
        T::Integer=200,
        prenormalized::Bool=false,
        signal_mask::Union{Nothing,AbstractVector}=nothing,
        padding_mask::Union{Nothing,AbstractVector}=nothing,
        history::Union{Nothing,Tuple}=nothing,
        h_mode_threshold::Union{Nothing,Real}=nothing)

    # ── Inputs ──
    hist, hmask, aux = history === nothing ? mean_normalized_history() : history
    size(hist) == (1, _MSE_HISTORY_LENGTH, _MAX_HISTORY_WIDTH) ||
        error("history_stats must be (1, $_MSE_HISTORY_LENGTH, $_MAX_HISTORY_WIDTH), got $(size(hist))")
    size(hmask) == (1, _MSE_HISTORY_LENGTH) ||
        error("history_masks must be (1, $_MSE_HISTORY_LENGTH), got $(size(hmask))")
    size(aux) == (1, _AUX_DIM) || error("aux_features must be (1, $_AUX_DIM), got $(size(aux))")

    if sequences === nothing
        sequences_z = zeros(Float32, T, _FPE_SIGNAL_DIM)
    else
        size(sequences, 2) == _FPE_SIGNAL_DIM ||
            error("sequences must have $_FPE_SIGNAL_DIM channels, got $(size(sequences,2))")
        T = size(sequences, 1)
        sequences_z = prenormalized ? Matrix{Float32}(sequences) :
                      normalized_signals(nn, sequences)
    end
    smask = signal_mask === nothing ? ones(Float32, _FPE_SIGNAL_DIM) : Float32.(signal_mask)
    length(smask) == _FPE_SIGNAL_DIM || error("signal_mask must have length $_FPE_SIGNAL_DIM")
    pmask = padding_mask === nothing ? ones(Float32, T) : Float32.(padding_mask)
    length(pmask) == T || error("padding_mask length ($(length(pmask))) must equal T ($T)")

    hist_f = Float32.(hist)
    hmask_f = Float32.(hmask)
    aux_f = Float32.(aux)
    seq3 = reshape(sequences_z, 1, T, _FPE_SIGNAL_DIM)
    smask2 = reshape(smask, 1, _FPE_SIGNAL_DIM)
    pmask2 = reshape(pmask, 1, T)

    # ── Defaults ──
    nan_T = fill(NaN32, T)
    raw = Dict{String,Vector{Float32}}()
    te_ped = copy(nan_T)
    ti_ped = copy(nan_T)
    t_rot_ped = copy(nan_T)
    edens_ped = copy(nan_T)
    edens_norm = copy(nan_T)
    hmode_logit_seq = copy(nan_T)
    hmode_prob_seq = copy(nan_T)
    is_h_mode_seq = falses(T)
    units = Dict{Symbol,String}()
    machine_embed = zeros(Float32, 1, 512)
    aux_embed = zeros(Float32, 1, 64)
    used_threshold = h_mode_threshold === nothing ? 0.5f0 : Float32(h_mode_threshold)

    # ── Run each bundle ──
    for (slug, b) in nn.bundles
        h_sliced = b.history_features == _MAX_HISTORY_WIDTH ? hist_f :
                   Float32.(@view hist_f[:, :, 1:b.history_features])

        mse_out = b.mse(Dict(
            "history_stats" => h_sliced,
            "history_masks" => hmask_f,
            "aux_features" => aux_f))
        me = mse_out["machine_embed"]::AbstractMatrix
        ae = mse_out["aux_embed"]::AbstractMatrix
        machine_embed = Matrix{Float32}(me)
        aux_embed = Matrix{Float32}(ae)

        fpe_out = b.fpe(Dict(
            "sequences" => seq3,
            "machine_embed" => Float32.(me),
            "aux_embed" => Float32.(ae),
            "signal_masks" => smask2,
            "padding_mask" => pmask2))
        pred_arr = fpe_out["predictions"]
        pred = _squeeze_predictions(pred_arr, T)
        raw[slug] = pred

        if is_regression(b)
            phys = @. pred * b.target_std + b.target_mean
            if b.target == "te_ped"
                te_ped = phys
                units[:te_ped] = b.target_units
            elseif b.target == "ti_ped"
                ti_ped = phys
                units[:ti_ped] = b.target_units
            elseif b.target == "t_rot_ped"
                t_rot_ped = phys
                units[:t_rot_ped] = b.target_units
            elseif b.target == "edens_ped"
                edens_ped = phys
                edens_norm = pred
                units[:edens_ped] = isempty(b.target_units) ? "1e19 m^-3" : b.target_units
            else
                @warn "PedestalNN: regression bundle $slug has unrecognised target=$(b.target); only available via raw[\"$slug\"]"
            end
        else
            thr = h_mode_threshold === nothing ? b.default_threshold : Float32(h_mode_threshold)
            used_threshold = thr
            hmode_logit_seq = pred
            hmode_prob_seq = Float32.(@. 1f0 / (1f0 + exp(-pred)))
            is_h_mode_seq = hmode_prob_seq .>= thr
        end
    end

    # Scalar summaries used by ActorPedestal (mean over the FPE window matches
    # the convention applied to predictions_physical further downstream).
    hmode_logit_scalar = any(isnan, hmode_logit_seq) ? NaN32 :
                         Float32(sum(hmode_logit_seq) / length(hmode_logit_seq))
    hmode_prob_scalar = any(isnan, hmode_prob_seq) ? NaN32 :
                        Float32(sum(hmode_prob_seq) / length(hmode_prob_seq))
    is_h_mode_scalar = isnan(hmode_prob_scalar) ? false : (hmode_prob_scalar >= used_threshold)

    return (;
        te_ped, ti_ped, t_rot_ped, edens_ped,
        hmode_logit_seq, hmode_prob_seq, is_h_mode_seq,
        hmode_logit=hmode_logit_scalar, hmode_prob=hmode_prob_scalar, is_h_mode=is_h_mode_scalar,
        predictions_physical=edens_ped, predictions_norm=edens_norm,
        raw, units, machine_embed, aux_embed)
end

"""
    predict_density(nn::PedestalNN; kwargs...) -> NamedTuple

Backward-compatible alias for [`predict_pedestal`](@ref). Returns the same
NamedTuple, including density (`predictions_physical`) along with the
temperature / rotation / H-mode fields.
"""
predict_density(nn::PedestalNN; kwargs...) = predict_pedestal(nn; kwargs...)

# Helper: ONNX FPE output may be (1, T) or (1, T, 1) depending on bundle export.
function _squeeze_predictions(arr, T::Integer)
    if ndims(arr) == 2
        size(arr, 1) == 1 && size(arr, 2) == T ||
            error("predictions tensor expected (1, $T), got $(size(arr))")
        return Float32.(vec(arr))
    elseif ndims(arr) == 3
        size(arr, 1) == 1 && size(arr, 2) == T && size(arr, 3) == 1 ||
            error("predictions tensor expected (1, $T, 1), got $(size(arr))")
        return Float32.(reshape(arr, T))
    else
        error("predictions tensor has unsupported ndims=$(ndims(arr)), shape=$(size(arr))")
    end
end
