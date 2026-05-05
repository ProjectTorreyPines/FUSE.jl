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
import NPZ
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

const _AUX_HISTORY_KEY = :nn_history_buffer

"""
    AUX_DEFAULT_SENTINEL :: NTuple{3,Float32}

Default `aux_features = [bzn_seconds, disrupt_seconds, disrupt_coverage]` triple
matching `HistoryManager._init_from_raw`'s sentinel branch (when no
`AuxFeaturesProvider` is wired). Per `pedestal-predictor-onnx/inference/aux_features.py`:

- `bzn_seconds = 0.0` — "no boronization data" sentinel.
- `disrupt_seconds = -1.0` — out-of-coverage sentinel (DIII-D shots > 180,847,
  i.e. the entire edensfit89 training shot range 196,000–205,000).
- `disrupt_coverage = 0.0` — false (correct for the same shot range).

Using `(0.0, -1.0, 0.0)` rather than `(0, 0, 0)` is closer to the H5/log path
the model was trained against, so mean-input runs don't accidentally tell the
encoder "this shot is in disrupt-list coverage with 0 s since last disruption".
"""
const AUX_DEFAULT_SENTINEL = (0f0, -1f0, 0f0)

"""
    mean_normalized_history() -> (hist, hmask, aux)

Produce the "neutral" MSE inputs corresponding to the normalization means
(z-score zeros for `history_stats`) plus the [`AUX_DEFAULT_SENTINEL`](@ref)
triple for `aux_features`. Use this when no machine-state history file is
available. Shapes match the v2 ONNX graph: `(1, 50, 458)`, `(1, 50)`, `(1, 3)`.
v1 bundles internally see `[..., 1:446]`.
"""
function mean_normalized_history()
    hist = zeros(Float32, 1, _MSE_HISTORY_LENGTH, _MAX_HISTORY_WIDTH)
    hmask = ones(Float32, 1, _MSE_HISTORY_LENGTH)
    aux = reshape(Float32[AUX_DEFAULT_SENTINEL...], 1, _AUX_DIM)
    return hist, hmask, aux
end

"""
    load_history_npz(path; norm_params_path=nothing) -> (hist, hmask, aux)

Load an MSE-history NPZ produced by PedestalPredictor's
`HistoryManager.save_state` (i.e. via `run_inference.py --save-history`),
returning a triple shaped to match the ONNX graph inputs and ready to pass
into [`predict_pedestal`](@ref) as the `history` keyword.

NPZ contract (per `pedestal-predictor-onnx/inference/history_manager.py::save_state`):

| key                  | dtype     | shape           | meaning                              |
|----------------------|-----------|-----------------|--------------------------------------|
| `history_stats`      | float32   | `(n_slots, F)`  | per-slot rolled-up signals (RAW)     |
| `history_masks`      | float32   | `(n_slots,)`    | 1.0 = real slot, 0.0 = padding       |
| `shot_ids`           | int64     | `(n_slots,)`    | shot id per slot (informational)     |
| `bzn_seconds`        | float     | `(1,)`          | seconds since last boronization      |
| `disrupt_seconds`    | float     | `(1,)`          | seconds since last disruption (-1 = sentinel) |
| `disrupt_coverage`   | float     | `(1,)`          | 1 if shot in disrupt-list coverage   |
| `n_history_features` | int       | `(1,)`          | 446 (v1) or 458 (v2)                 |

Returned arrays:
- `hist::Array{Float32,3}`  — `(1, _MSE_HISTORY_LENGTH, _MAX_HISTORY_WIDTH)` z-scored
   if `norm_params_path` was provided, else raw.
- `hmask::Matrix{Float32}`  — `(1, _MSE_HISTORY_LENGTH)`
- `aux::Matrix{Float32}`    — `(1, 3)` = `[bzn_sec, dis_sec, dis_cov]` (raw, not z-scored)

Adaptation behavior:
- Slot count != `_MSE_HISTORY_LENGTH (=50)` → zero-pad (stats and mask) or take
  the **latest** `_MSE_HISTORY_LENGTH` slots.
- Width < `_MAX_HISTORY_WIDTH (=458)` (v1_446 NPZs) → zero-pad on the right; v1
  bundles ignore the trailing 12 features via the slicing in `predict_pedestal`.
- `norm_params_path` should point to a JSON with either
  `history_means`/`history_stds` or `means`/`stds` keys (length F). Without it,
  raw stats are returned and a warning is emitted — useful only for offline
  debugging since the encoder was trained on z-scored inputs.

Mirrors what `HistoryManager.get_mse_inputs` produces in Python so this is the
on-ramp for `--save-history` round-trip / parity testing.
"""
function load_history_npz(path::AbstractString;
        norm_params_path::Union{Nothing,AbstractString}=nothing)
    isfile(path) || error("load_history_npz: NPZ not found at $path")
    data = NPZ.npzread(path)
    return _assemble_history(data; norm_params_path, source="NPZ at $path")
end

"""
    _assemble_history(record::AbstractDict; norm_params_path=nothing, source="<source>") -> (hist, hmask, aux)

Shared back-end for both [`load_history_npz`](@ref) and
[`mse_history_from_aux`](@ref). Validates the required keys, optionally
z-scores `history_stats`, pads/truncates to the encoder's `(1, 50, 458)`
contract, and assembles the `(1, 3)` aux triple from `[bzn, dis, cov]`.

`record` accepts any AbstractDict (or NamedTuple-like) with the same key set
as PedestalPredictor's `HistoryManager.save_state` NPZ. `source` is used in
error messages only.
"""
function _assemble_history(record;
        norm_params_path::Union{Nothing,AbstractString}=nothing,
        source::AbstractString="<unknown>")
    # NamedTuple only supports Symbol keys; AbstractDict from NPZ uses String.
    # Try Symbol first (cheap, total for NamedTuple), then fall back to String.
    function _get(k)
        sym = Symbol(k)
        if record isa NamedTuple
            return haskey(record, sym) ? getfield(record, sym) : nothing
        end
        return haskey(record, k) ? record[k] :
               (haskey(record, sym) ? record[sym] : nothing)
    end
    for k in ("history_stats", "history_masks", "bzn_seconds",
              "disrupt_seconds", "disrupt_coverage")
        _get(k) === nothing && error("_assemble_history: $source missing required key '$k'")
    end

    raw_stats = _get("history_stats")
    raw_masks = _get("history_masks")
    ndims(raw_stats) == 2 || error("_assemble_history: history_stats must be 2-D, got $(ndims(raw_stats))-D")
    ndims(raw_masks) == 1 || error("_assemble_history: history_masks must be 1-D, got $(ndims(raw_masks))-D")
    size(raw_stats, 1) == size(raw_masks, 1) ||
        error("_assemble_history: history_stats rows ($(size(raw_stats,1))) != history_masks length ($(size(raw_masks,1)))")

    n_slots_in = size(raw_stats, 1)
    n_features = size(raw_stats, 2)
    n_features <= _MAX_HISTORY_WIDTH ||
        error("_assemble_history: history width $n_features exceeds max $_MAX_HISTORY_WIDTH")

    n_hist_decl_v = _get("n_history_features")
    if n_hist_decl_v !== nothing
        n_hist_decl = Int(first(n_hist_decl_v))
        n_hist_decl == n_features ||
            error("_assemble_history: n_history_features=$n_hist_decl != actual width $n_features")
        n_hist_decl in (446, 458) ||
            @warn "_assemble_history: unusual history width $n_hist_decl (expected 446 or 458)"
    end

    stats_typed = Float32.(raw_stats)
    if norm_params_path !== nothing
        isfile(norm_params_path) ||
            error("_assemble_history: norm_params_path not found: $norm_params_path")
        params = JSON.parsefile(norm_params_path)
        # Match Python `HistoryManager._load_norm_params`: only the per-feature
        # `history_means`/`history_stds` keys (length == 446 or 458) trigger
        # z-scoring. Per-base-signal files (e.g. `normalization_params_prmtan67.json`,
        # which has `means`/`stds` of length 67) are intentionally NOT accepted
        # here — they encode a different normalization step that's already
        # baked into `compute_shot_stats`'s raw-signal preprocessing.
        hmean = get(params, "history_means", get(params, "history_mean", nothing))
        hstd  = get(params, "history_stds",  get(params, "history_std",  nothing))
        if hmean === nothing || hstd === nothing
            @warn "_assemble_history: $norm_params_path lacks history_means/history_stds — history will be passed un-normalized"
        else
            length(hmean) == n_features ||
                error("_assemble_history: history_means length $(length(hmean)) != $n_features")
            length(hstd) == n_features ||
                error("_assemble_history: history_stds length $(length(hstd)) != $n_features")
            μ = reshape(Float32.(hmean), 1, n_features)
            σ = reshape(Float32.(hstd),  1, n_features)
            @. stats_typed = (stats_typed - μ) / (σ + 1f-8)
        end
    else
        @warn "_assemble_history: no norm_params_path supplied for $source — history will be passed un-normalized (matches Python `HistoryManager` default when no normalization_params.json contains history_means/history_stds)"
    end

    stats_padded = zeros(Float32, _MSE_HISTORY_LENGTH, _MAX_HISTORY_WIDTH)
    masks_padded = zeros(Float32, _MSE_HISTORY_LENGTH)
    n_take = min(n_slots_in, _MSE_HISTORY_LENGTH)
    src_lo = n_slots_in - n_take + 1
    dst_lo = _MSE_HISTORY_LENGTH - n_take + 1
    stats_padded[dst_lo:end, 1:n_features] .= @view stats_typed[src_lo:end, :]
    masks_padded[dst_lo:end] .= @view Float32.(raw_masks)[src_lo:end]

    hist = reshape(stats_padded, 1, _MSE_HISTORY_LENGTH, _MAX_HISTORY_WIDTH)
    hmask = reshape(masks_padded, 1, _MSE_HISTORY_LENGTH)

    aux = reshape(Float32[
        Float32(first(_get("bzn_seconds"))),
        Float32(first(_get("disrupt_seconds"))),
        Float32(first(_get("disrupt_coverage"))),
    ], 1, _AUX_DIM)

    return hist, hmask, aux
end

# ── Live history buffer (dd._aux[:nn_history_buffer]) ───────────────────────
#
# Schema (NamedTuple identical to HistoryManager.save_state's NPZ):
#   history_stats        :: Matrix{Float32}   # (n_slots, n_features) — RAW
#   history_masks        :: Vector{Float32}   # (n_slots,)
#   shot_ids             :: Vector{Int64}     # (n_slots,)
#   bzn_seconds          :: Vector{Float64}   # (1,)
#   disrupt_seconds      :: Vector{Float64}   # (1,)
#   disrupt_coverage     :: Vector{Float64}   # (1,)
#   n_history_features   :: Vector{Int64}     # (1,) — 446 or 458
#
# The matching schema means save_history_npz / load_history_npz round-trip
# directly with this struct, and `mse_history_from_aux` uses the same
# `_assemble_history` back-end as `load_history_npz`.

"""
    init_history_buffer(; n_features=458, bzn_seconds=0.0, disrupt_seconds=-1.0,
                        disrupt_coverage=0.0) -> NamedTuple

Construct an empty live MSE-history buffer matching `HistoryManager.save_state`'s
NPZ schema. Use `n_features = 446` for v1-only deployments. Aux scalars default
to [`AUX_DEFAULT_SENTINEL`](@ref).
"""
function init_history_buffer(;
        n_features::Integer=_MAX_HISTORY_WIDTH,
        bzn_seconds::Real=AUX_DEFAULT_SENTINEL[1],
        disrupt_seconds::Real=AUX_DEFAULT_SENTINEL[2],
        disrupt_coverage::Real=AUX_DEFAULT_SENTINEL[3])
    n_features in (446, 458) ||
        @warn "init_history_buffer: unusual n_features=$n_features (expected 446 or 458)"
    return (
        history_stats      = Matrix{Float32}(undef, 0, Int(n_features)),
        history_masks      = Float32[],
        shot_ids           = Int64[],
        bzn_seconds        = Float64[bzn_seconds],
        disrupt_seconds    = Float64[disrupt_seconds],
        disrupt_coverage   = Float64[disrupt_coverage],
        n_history_features = Int64[Int(n_features)],
    )
end

"""
    push_shot_history!(dd, stats_vec; shot_id, mask=1.0, max_slots=50,
                       bzn_seconds=nothing, disrupt_seconds=nothing,
                       disrupt_coverage=nothing) -> NamedTuple

Append a shot-summary stats row to `dd._aux[:nn_history_buffer]`. Creates the
buffer if absent (using `length(stats_vec)` as the feature width). The buffer
is kept bounded at `max_slots` rows by dropping the oldest entries. Aux
scalars are updated in place when supplied; otherwise their previous values
(or sentinel defaults) are preserved.

The expected workflow at end-of-shot:

```julia
stats = compute_shot_stats(dd)         # Vector{Float32}, length 446 or 458 — TODO
push_shot_history!(dd, stats; shot_id=current_shot)
```

so the next shot's `predict_pedestal` call sees this shot's contribution via
[`mse_history_from_aux`](@ref).
"""
function push_shot_history!(dd, stats_vec::AbstractVector;
        shot_id::Integer,
        mask::Real=1f0,
        max_slots::Integer=_MSE_HISTORY_LENGTH,
        bzn_seconds::Union{Nothing,Real}=nothing,
        disrupt_seconds::Union{Nothing,Real}=nothing,
        disrupt_coverage::Union{Nothing,Real}=nothing)
    aux = getfield(dd, :_aux)
    if !haskey(aux, _AUX_HISTORY_KEY)
        aux[_AUX_HISTORY_KEY] = init_history_buffer(; n_features=length(stats_vec))
    end
    buf = aux[_AUX_HISTORY_KEY]
    n_features = Int(first(buf.n_history_features))
    length(stats_vec) == n_features ||
        error("push_shot_history!: stats_vec length $(length(stats_vec)) != buffer width $n_features")

    new_row = reshape(Float32.(stats_vec), 1, n_features)
    new_stats = vcat(buf.history_stats, new_row)
    new_masks = vcat(buf.history_masks, Float32(mask))
    new_ids   = vcat(buf.shot_ids, Int64(shot_id))

    # Drop oldest entries if over capacity
    if size(new_stats, 1) > max_slots
        cut = size(new_stats, 1) - Int(max_slots)
        new_stats = new_stats[cut+1:end, :]
        new_masks = new_masks[cut+1:end]
        new_ids   = new_ids[cut+1:end]
    end

    bzn_v = bzn_seconds === nothing ? buf.bzn_seconds : Float64[bzn_seconds]
    dis_v = disrupt_seconds === nothing ? buf.disrupt_seconds : Float64[disrupt_seconds]
    cov_v = disrupt_coverage === nothing ? buf.disrupt_coverage : Float64[disrupt_coverage]

    aux[_AUX_HISTORY_KEY] = (
        history_stats      = new_stats,
        history_masks      = new_masks,
        shot_ids           = new_ids,
        bzn_seconds        = bzn_v,
        disrupt_seconds    = dis_v,
        disrupt_coverage   = cov_v,
        n_history_features = buf.n_history_features,
    )
    return aux[_AUX_HISTORY_KEY]
end

"""
    mse_history_from_aux(dd, nn::PedestalNN; norm_params_path=nothing) -> (hist, hmask, aux)

Build the `(hist, hmask, aux)` triple that `predict_pedestal` consumes from
the MSE history buffer stored in `dd._aux[:nn_history_buffer]`. Returns
[`mean_normalized_history`](@ref) when the buffer is absent or empty, so this
is a drop-in replacement for the trained-mean fallback.

# Data-source contract

The MSE history is **not** supplied by GSLite on the real-time ZMQ wire —
GSLite only provides the 32-channel FPE actuator state via `dd._aux[:zmq_*]`.
The 50-slot × 446/458-feature history tensor captures a slow machine-state
context (boronization age, last-week coil usage, gas mix, etc.) that
requires fetching ~67 PTDATA signals from MDSplus per shot and rolling
them up offline. It reaches `dd._aux[:nn_history_buffer]` via one of:

1. **Static NPZ** (production path today) — loaded once at FUSE startup
   from `ENV["FUSE_PEDESTAL_NN_HISTORY_NPZ"]` (e.g. the reference snapshot
   in `FUSE/data/pedestal_nn/history_morning.npz`) via
   [`load_history_npz`](@ref). Operator-refreshed between run-days.
2. **Live FUSE-side bookkeeping** (future) — `compute_shot_stats(dd)` +
   [`push_shot_history!`](@ref) at the end of each shot in FUSE itself.
   Currently stubbed; wires up the feed-forward path so this function
   becomes self-sufficient once the per-block source registries land.

`norm_params_path` should point at the encoder's history-norm JSON (typically
`onnx_models/normalization_params_prmtan67.json` or its `prmtan_clean`
variant). Without it, raw stats pass through unscaled and a warning is emitted
— useful for debugging but the encoder was trained on z-scored inputs.
"""
function mse_history_from_aux(dd, nn::PedestalNN;
        norm_params_path::Union{Nothing,AbstractString}=nothing)
    aux = getfield(dd, :_aux)
    if !haskey(aux, _AUX_HISTORY_KEY) || isempty(getproperty(aux[_AUX_HISTORY_KEY], :history_masks))
        return mean_normalized_history()
    end
    return _assemble_history(aux[_AUX_HISTORY_KEY]; norm_params_path,
                             source="dd._aux[:$(_AUX_HISTORY_KEY)]")
end

"""
    save_history_npz(dd, path) -> path

Dump `dd._aux[:nn_history_buffer]` to `path` as an NPZ matching
`HistoryManager.save_state`'s schema (round-trippable with
[`load_history_npz`](@ref)). Errors if the buffer is absent.
"""
function save_history_npz(dd, path::AbstractString)
    aux = getfield(dd, :_aux)
    haskey(aux, _AUX_HISTORY_KEY) ||
        error("save_history_npz: dd._aux[:$(_AUX_HISTORY_KEY)] missing; nothing to save")
    buf = aux[_AUX_HISTORY_KEY]
    NPZ.npzwrite(path, Dict(
        "history_stats"      => buf.history_stats,
        "history_masks"      => buf.history_masks,
        "shot_ids"           => buf.shot_ids,
        "bzn_seconds"        => buf.bzn_seconds,
        "disrupt_seconds"    => buf.disrupt_seconds,
        "disrupt_coverage"   => buf.disrupt_coverage,
        "n_history_features" => buf.n_history_features,
    ))
    return path
end

"""
    compute_shot_stats(dd; dataset_version="v2_458") -> Vector{Float32}

Compute the per-shot stats vector (length 446 for v1_446, 458 for v2_458) for
the current shot in `dd`, suitable for [`push_shot_history!`](@ref). Mirrors
the Python `HistoryManager.compute_shot_stats` layout (Block 1–6 from
`pedestal-predictor-onnx/inference/history_manager.py` and
`docs/onnx_export_spec.md` §"Block 1: Base signals (indices 0-267, 67 signals
x 4)").

!!! warning "Stub"
    This is a skeleton: the per-block source registries that map each of the
    67 Block-1 signals (kappa, betap, …, density, n1rms, n2rms, prad), the 40
    ECE channels (Block 3), the radiation diagnostics (Block 4), and the
    self-referential blocks 2/5/6 (previous-shot edensfit100/edensfit89/te_ped/
    ti_ped/t_rot_ped) onto FUSE `dd` paths and `dd._aux[:zmq_*]` records are
    **not yet wired**. This stub returns all-zero stats with a warning so the
    buffer plumbing can be exercised end-to-end. Filling in the registries is
    the next milestone after this skeleton lands.
"""
function compute_shot_stats(dd; dataset_version::AbstractString="v2_458")
    haskey(_DATASET_WIDTH, dataset_version) ||
        error("compute_shot_stats: unknown dataset_version=$dataset_version (expected v1_446 or v2_458)")
    n = _DATASET_WIDTH[dataset_version]
    @warn "compute_shot_stats: stub returning all-zero stats — Block 1–6 source registries not yet wired"
    return zeros(Float32, n)
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
