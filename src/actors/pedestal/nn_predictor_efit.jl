#= ============================================================================ =#
#  fuse_efit — DIII-D Mamba-2 pedestal predictor with EFIT geometry inputs (ONNX) #
#                                                                                #
#  Upstream name FUSE-EFIT. Published at https://huggingface.co/SCS-Lab/FUSE-PED #
#  Source repo: https://github.com/NoahH72/FUSE-Pedestal-Predictor-Mamba-ONNX    #
#  (branch `fuse-efit`; docs/EXPORT.md, docs/INPUT_CONTRACT.md, HANDOFF.md)      #
#                                                                                #
#  Successor of fuse29 (nn_predictor.jl) with a different input set: the 20      #
#  coil currents are gone; the measured plasma current, 14 boundary-shape        #
#  quantities of the equilibrium and the loop voltage come in.                   #
#                                                                                #
#    inputs      (batch, 25)  physical units, FUSE_EFIT_INPUT_NAMES order        #
#    indicators  (batch, 5)   0/1 per tick: efit_valid, vloop_valid,             #
#                             rvsout_present, rvsou_present, rvsiu_present       #
#    available   (batch, 25)  0/1 per SHOT and channel (constant within a shot)  #
#    conv_state  (L, batch, 1792, 3), ssm_state (L, batch, 12, 128, 64)          #
#                             L = 6 (form :single) or 48 (form :averaged)        #
#    predictions (batch, 9)   same nine heads as fuse29, physical units; heads   #
#                             5 and 6 are called ne_knee / te_knee here          #
#                                                                                #
#  Normalisation, the ±20σ clip and the indicator/availability masking are       #
#  inside the graphs: a column whose governing indicator or availability flag    #
#  is 0 may hold anything (NaN included). NaN in a column that IS in use is      #
#  silently read as the training mean, so never pass NaN for a value you have.   #
#                                                                                #
#  Two forms with the same calling convention: `:single` (one network, the       #
#  frozen delivered model, ~5 ms per tick) and `:averaged` (eight networks       #
#  averaged at every tick: lower error, ~10x the cost, 1.5 GB of graphs).        #
#                                                                                #
#  Artefacts live *outside* FUSE. Their root (the directory holding              #
#  `manifest.json` and `era/<form>/`) is resolved from, in order:                #
#    (1) the `onnx_dir` kwarg passed to `load_pedestal_efit_nn`,                 #
#    (2) the environment variable `FUSE_PEDESTAL_EFIT_NN_DIR`,                   #
#    (3) `<pedestal-predictor-mamba-onnx-fuse repo>/artifacts/fuse_efit/`.       #
#  Missing bundles are auto-fetched from HuggingFace (the requested form only).  #
#= ============================================================================ =#

const PEDESTAL_EFIT_NN_ENV = "FUSE_PEDESTAL_EFIT_NN_DIR"
const PEDESTAL_EFIT_NN_HF_REPO_ENV = "FUSE_PEDESTAL_EFIT_NN_HF_REPO"
const PEDESTAL_EFIT_NN_HF_REVISION_ENV = "FUSE_PEDESTAL_EFIT_NN_HF_REVISION"
const PEDESTAL_EFIT_NN_HF_REPO_DEFAULT = "SCS-Lab/FUSE-PED"
const PEDESTAL_EFIT_NN_HF_REVISION_DEFAULT = "main"
const PEDESTAL_EFIT_NN_PROTOCOL = "era"   # the delivered (in-era) model; the bundle is `era/<form>/`

# ── The fuse_efit I/O contract (mirrors model_config.json; asserted against the bundle at load) ──
const FUSE_EFIT_INPUT_NAMES = (
    "pinj", "tinj", "ech_mw",
    "gasa_cal", "gasb_cal", "gasc_cal", "gasd_cal", "gase_cal",
    "bt", "ip",
    "aminor", "rsurf", "zsurf", "kappa", "tritop", "tribot", "drsep",
    "rvsout", "gapin", "rvsou", "rvsiu", "shafranov", "zeta_ou", "zeta_ol",
    "vloop")
const FUSE_EFIT_INPUT_UNITS = (
    "MW", "N*m", "MW",
    "Torr*L/s", "Torr*L/s", "Torr*L/s", "Torr*L/s", "Torr*L/s",
    "T", "MA",
    "m", "m", "m", "-", "-", "-", "m",
    "m", "m", "m", "m", "m", "-", "-",
    "V")
const FUSE_EFIT_INDICATOR_NAMES = ("efit_valid", "vloop_valid", "rvsout_present", "rvsou_present", "rvsiu_present")
# the 14 boundary-shape channels governed by `efit_valid`
const FUSE_EFIT_SHAPE_NAMES = ("aminor", "rsurf", "zsurf", "kappa", "tritop", "tribot", "drsep",
                               "rvsout", "gapin", "rvsou", "rvsiu", "shafranov", "zeta_ou", "zeta_ol")
const FUSE_EFIT_HEAD_NAMES = ("ne", "te_ped", "ti_ped", "t_rot_ped",
                              "ne_knee", "te_knee", "rho_sym", "ne_top_loc",
                              "hmode")
const FUSE_EFIT_HEAD_UNITS = ("1e19 m^-3", "keV", "keV", "krad/s",
                              "1e19 m^-3", "keV", "rho_tor", "rho_tor",
                              "probability")
const FUSE_EFIT_GATED_HEADS = ("ne_knee", "te_knee", "rho_sym", "ne_top_loc")
const FUSE_EFIT_N_INPUTS = length(FUSE_EFIT_INPUT_NAMES)
const FUSE_EFIT_N_INDICATORS = length(FUSE_EFIT_INDICATOR_NAMES)
const FUSE_EFIT_N_HEADS = length(FUSE_EFIT_HEAD_NAMES)
const FUSE_EFIT_STEP_INPUTS = ("inputs", "indicators", "available", "conv_state", "ssm_state")
const FUSE_EFIT_STEP_OUTPUTS = ("predictions", "conv_state_out", "ssm_state_out")
const FUSE_EFIT_SEQ_INPUTS = ("inputs", "indicators", "available")
const FUSE_EFIT_FORMS = (:single, :averaged)
const _FUSE_EFIT_TOPLEVEL_FILES = ("manifest.json", "SHA256SUMS")
const FUSE_EFIT_RHO_REF = 0.85   # radius the four continuous heads are defined at (model_config outputs.predictions.meaning)

"""
    fuse_efit_bundle_files(form::Symbol) -> NTuple{6,String}

File names of one fuse_efit bundle (`era/<form>/`).
"""
function fuse_efit_bundle_files(form::Symbol)
    form in FUSE_EFIT_FORMS || error("fuse_efit: form must be one of $(FUSE_EFIT_FORMS), got $(repr(form))")
    return ("fuse_efit_$(form)_step.onnx", "fuse_efit_$(form)_sequence.onnx", "model_config.json",
            "norms.json", "input_ranges.json", "provenance.json")
end

"""
    FuseEfitNN

One loaded fuse_efit bundle (one form): the stateful step graph, the fixed-length
sequence graph (loaded on first use) and the JSON sidecars.

Fields:
- `step`               : ONNX Runtime session of the step graph
- `seq`                : `Ref` to the sequence-graph session (`nothing` until [`fuse_efit_sequence`](@ref) needs it)
- `config`             : parsed `model_config.json`
- `input_names`        : 25 input channel names, in graph column order
- `indicator_names`    : 5 per-tick indicator names, in graph column order
- `head_names`         : 9 output head names, in graph column order
- `input_ranges`       : parsed `input_ranges.json` (training-split distribution, or `nothing`)
- `seq_len`            : fixed sequence length of the sequence graph (256)
- `form`               : `:single` or `:averaged`
- `hmode_balanced`, `hmode_conservative` : the form's two validation-set gate thresholds
- `period_s`           : model cadence (0.05 s)
- `rho_ref`            : radius the continuous heads are defined at (rho_tor = 0.85)
- `conv_shape`, `ssm_shape` : recurrent state shapes for batch 1
- `bundle_dir`         : directory the artefacts were loaded from
"""
struct FuseEfitNN
    step::ONNXRunTime.InferenceSession
    seq::Base.RefValue{Union{Nothing,ONNXRunTime.InferenceSession}}
    config::Dict{String,Any}
    input_names::Vector{String}
    indicator_names::Vector{String}
    head_names::Vector{String}
    input_ranges::Union{Nothing,Dict{String,Any}}
    seq_len::Int
    form::Symbol
    hmode_balanced::Float64
    hmode_conservative::Float64
    period_s::Float64
    rho_ref::Float64
    conv_shape::NTuple{4,Int}
    ssm_shape::NTuple{5,Int}
    bundle_dir::String
end

function Base.show(io::IO, nn::FuseEfitNN)
    print(io, "FuseEfitNN(form=$(nn.form), $(FUSE_EFIT_N_INPUTS)+$(FUSE_EFIT_N_INDICATORS)→$(FUSE_EFIT_N_HEADS), period=$(nn.period_s)s, rho_ref=$(nn.rho_ref), dir=\"$(nn.bundle_dir)\")")
end

# ── Path resolution ─────────────────────────────────────────────────────────

"""
    resolve_pedestal_efit_nn_dir(; onnx_dir=nothing) -> String

Resolve the fuse_efit artefact root: the directory that holds `manifest.json`
and `era/single/`, `era/averaged/` (the layout of the HuggingFace download).
`load_pedestal_efit_nn` also accepts the `era/` directory itself or a bare
bundle directory containing `model_config.json`.
"""
function resolve_pedestal_efit_nn_dir(; onnx_dir::Union{Nothing,AbstractString}=nothing)
    if onnx_dir !== nothing
        return String(onnx_dir)
    end
    if haskey(ENV, PEDESTAL_EFIT_NN_ENV)
        return ENV[PEDESTAL_EFIT_NN_ENV]
    end
    return normpath(joinpath(__FUSE__, "..", "pedestal-predictor-mamba-onnx-fuse", "artifacts", "fuse_efit"))
end

function _fuse_efit_bundle_complete(bdir::AbstractString, form::Symbol)
    isdir(bdir) || return false
    for f in fuse_efit_bundle_files(form)
        path = joinpath(bdir, f)
        (isfile(path) && filesize(path) > 0) || return false
    end
    return true
end

# A bare bundle dir wins, then `<root>/<form>/` (root is the `era` folder), then `<root>/era/<form>/`.
function _fuse_efit_find_bundle_dir(root::AbstractString, form::Symbol)
    for cand in (String(root), joinpath(root, String(form)), joinpath(root, PEDESTAL_EFIT_NN_PROTOCOL, String(form)))
        _fuse_efit_bundle_complete(cand, form) && return cand
    end
    return nothing
end

# ── HuggingFace auto-fetch ──────────────────────────────────────────────────

"""
    pedestal_efit_nn_files(; form=:single) -> Vector{String}

Relative paths (under the artefact root) that [`download_pedestal_efit_nn!`](@ref)
fetches from HuggingFace for one form: `manifest.json`, `SHA256SUMS` and the six
files of `era/<form>/`.
"""
function pedestal_efit_nn_files(; form::Symbol=:single)
    files = String[String(f) for f in _FUSE_EFIT_TOPLEVEL_FILES]
    for f in fuse_efit_bundle_files(form)
        push!(files, "$(PEDESTAL_EFIT_NN_PROTOCOL)/$(form)/$f")
    end
    return files
end

"""
    download_pedestal_efit_nn!(target_dir::AbstractString;
                               repo=PEDESTAL_EFIT_NN_HF_REPO_DEFAULT,
                               revision=PEDESTAL_EFIT_NN_HF_REVISION_DEFAULT,
                               form=:single, force=false, verbose=true) -> target_dir

Pull one fuse_efit form from HuggingFace into `target_dir`, producing the layout
`load_pedestal_efit_nn` expects (`manifest.json` + `era/<form>/…`). Idempotent:
files already present are skipped unless `force=true`. Every file's size is
checked against `manifest.json` after the download (a truncated graph would
otherwise fail much later, inside ONNX Runtime).

`:single` is ~184 MB (two ~92 MB graphs), `:averaged` ~1.5 GB. `repo` /
`revision` are also overridable via `ENV[$(repr(PEDESTAL_EFIT_NN_HF_REPO_ENV))]`
and `ENV[$(repr(PEDESTAL_EFIT_NN_HF_REVISION_ENV))]`. `load_pedestal_efit_nn`
calls this automatically when the requested form is missing; opt out with
`ENV["$PEDESTAL_NN_AUTOFETCH_ENV"] = "0"` (shared with fuse29).
"""
function download_pedestal_efit_nn!(target_dir::AbstractString;
        repo::AbstractString=get(ENV, PEDESTAL_EFIT_NN_HF_REPO_ENV, PEDESTAL_EFIT_NN_HF_REPO_DEFAULT),
        revision::AbstractString=get(ENV, PEDESTAL_EFIT_NN_HF_REVISION_ENV, PEDESTAL_EFIT_NN_HF_REVISION_DEFAULT),
        form::Symbol=:single,
        force::Bool=false,
        verbose::Bool=true)
    mkpath(target_dir)
    files = pedestal_efit_nn_files(; form)
    base_url = "https://huggingface.co/$repo/resolve/$revision"
    verbose && @info "fuse_efit: fetching $(length(files)) files (form $form) from $repo @ $revision -> $target_dir"
    for (i, rel) in enumerate(files)
        dest = joinpath(target_dir, rel)
        if !force && isfile(dest) && filesize(dest) > 0
            verbose && @info "  [$i/$(length(files))] cached $rel ($(filesize(dest)) B)"
            continue
        end
        url = "$base_url/$rel"
        mkpath(dirname(dest))
        verbose && @info "  [$i/$(length(files))] fetching $rel"
        try
            Downloads.download(url, dest)
        catch err
            isfile(dest) && rm(dest; force=true)  # clean up partial file
            error("download_pedestal_efit_nn!: failed to fetch $url\n  -> $err")
        end
        if filesize(dest) == 0
            rm(dest; force=true)
            error("download_pedestal_efit_nn!: $rel downloaded empty from $url")
        end
        verbose && @info "  [$i/$(length(files))] -> $(filesize(dest)) B"
    end
    for msg in fuse_efit_verify_sizes(target_dir; form)
        error("download_pedestal_efit_nn!: $msg (delete the file and retry, or pass force=true)")
    end
    verbose && @info "fuse_efit: download complete ($(length(files)) files)"
    return String(target_dir)
end

"""
    fuse_efit_verify_sizes(root::AbstractString; form=:single) -> Vector{String}

Compare the byte size of every file of `era/<form>/` under `root` with
`root/manifest.json`. Returns one message per mismatch; empty means the bundle
is complete (or there is no manifest to check against).
"""
function fuse_efit_verify_sizes(root::AbstractString; form::Symbol=:single)
    msgs = String[]
    mpath = joinpath(root, "manifest.json")
    isfile(mpath) || return msgs
    manifest = JSON.parsefile(mpath; dicttype=Dict{String,Any})
    listed = get(manifest, "files", Dict{String,Any}())
    for f in fuse_efit_bundle_files(form)
        rel = "$(PEDESTAL_EFIT_NN_PROTOCOL)/$(form)/$f"
        haskey(listed, rel) || continue
        path = joinpath(root, rel)
        want = Int(listed[rel]["bytes"])
        have = isfile(path) ? filesize(path) : 0
        have == want || push!(msgs, "$rel is $have bytes, manifest.json says $want")
    end
    return msgs
end

function _fuse_efit_maybe_autofetch!(root::AbstractString, form::Symbol)
    autofetch = lowercase(strip(get(ENV, PEDESTAL_NN_AUTOFETCH_ENV, "")))
    autofetch in ("0", "false", "no", "off") && return false
    @info "fuse_efit: ONNX bundle for form $form missing or incomplete; auto-fetching from HuggingFace ($(get(ENV, PEDESTAL_EFIT_NN_HF_REPO_ENV, PEDESTAL_EFIT_NN_HF_REPO_DEFAULT)) @ $(get(ENV, PEDESTAL_EFIT_NN_HF_REVISION_ENV, PEDESTAL_EFIT_NN_HF_REVISION_DEFAULT))) -> $root\n  Set ENV[\"$PEDESTAL_NN_AUTOFETCH_ENV\"] = \"0\" to opt out."
    download_pedestal_efit_nn!(root; form)
    return true
end

# ── Loading ─────────────────────────────────────────────────────────────────

const _PEDESTAL_EFIT_NN_CACHE = Dict{Tuple{String,Symbol},FuseEfitNN}()

"""
    load_pedestal_efit_nn(; onnx_dir=nothing, form=:single) -> FuseEfitNN

Load one fuse_efit bundle. `onnx_dir` (or `ENV["FUSE_PEDESTAL_EFIT_NN_DIR"]`)
may be the download root holding `manifest.json` and `era/`, the `era/`
directory, or a bare bundle directory containing `model_config.json`.
Cached per process; repeated calls return the same object.

The bundle is validated against the contract in this file at load time (input,
indicator and head names/order, graph tensor names, state shapes), so a
mismatched bundle fails here rather than producing plausible nonsense downstream.
"""
function load_pedestal_efit_nn(; onnx_dir::Union{Nothing,AbstractString}=nothing, form::Symbol=:single)
    form in FUSE_EFIT_FORMS || error("fuse_efit: form must be one of $(FUSE_EFIT_FORMS), got $(repr(form))")
    root = abspath(resolve_pedestal_efit_nn_dir(; onnx_dir))
    cache_key = (root, form)
    haskey(_PEDESTAL_EFIT_NN_CACHE, cache_key) && return _PEDESTAL_EFIT_NN_CACHE[cache_key]
    bdir = _fuse_efit_find_bundle_dir(root, form)
    if bdir === nothing
        _fuse_efit_maybe_autofetch!(root, form)
        bdir = _fuse_efit_find_bundle_dir(root, form)
    end
    bdir === nothing && error("fuse_efit: no complete `$form` bundle under $root\n" *
        "Expected $(collect(fuse_efit_bundle_files(form))) in `$root/$(PEDESTAL_EFIT_NN_PROTOCOL)/$form/`. Either:\n" *
        "  • Unset / set ENV[\"$PEDESTAL_NN_AUTOFETCH_ENV\"]=\"1\" and retry (downloads from HuggingFace)\n" *
        "  • Run `FUSE.download_pedestal_efit_nn!(\"$root\"; form=$(repr(form)))` explicitly\n" *
        "  • Set ENV[\"$PEDESTAL_EFIT_NN_ENV\"] or pass `onnx_dir=...` to point at an existing copy")
    nn = _load_fuse_efit(bdir, form)
    _PEDESTAL_EFIT_NN_CACHE[cache_key] = nn
    return nn
end

function _load_fuse_efit(bdir::AbstractString, form::Symbol)
    config = JSON.parsefile(joinpath(bdir, "model_config.json"); dicttype=Dict{String,Any})

    # -- contract assertions (names & order are load-bearing) --
    cfg_form = Symbol(get(config, "form", String(form)))
    cfg_form == form || error("fuse_efit: bundle in $bdir is form `$cfg_form` but `$form` was requested")
    input_names = String.(config["inputs"]["inputs"]["names"])
    indicator_names = String.(config["inputs"]["indicators"]["names"])
    available_names = String.(config["inputs"]["available"]["names"])
    head_names = String.(config["outputs"]["predictions"]["names"])
    gated = String.(config["outputs"]["hmode_gated_heads"])
    input_names == collect(FUSE_EFIT_INPUT_NAMES) ||
        error("fuse_efit: input names/order in $(bdir)/model_config.json do not match FUSE_EFIT_INPUT_NAMES:\n  bundle: $input_names\n  FUSE:   $(collect(FUSE_EFIT_INPUT_NAMES))")
    available_names == input_names ||
        error("fuse_efit: `available` names in $(bdir)/model_config.json differ from the input names: $available_names")
    indicator_names == collect(FUSE_EFIT_INDICATOR_NAMES) ||
        error("fuse_efit: indicator names/order in $(bdir)/model_config.json do not match FUSE_EFIT_INDICATOR_NAMES:\n  bundle: $indicator_names\n  FUSE:   $(collect(FUSE_EFIT_INDICATOR_NAMES))")
    head_names == collect(FUSE_EFIT_HEAD_NAMES) ||
        error("fuse_efit: head names/order in $(bdir)/model_config.json do not match FUSE_EFIT_HEAD_NAMES:\n  bundle: $head_names\n  FUSE:   $(collect(FUSE_EFIT_HEAD_NAMES))")
    gated == collect(FUSE_EFIT_GATED_HEADS) ||
        error("fuse_efit: hmode-gated heads in bundle ($gated) differ from FUSE_EFIT_GATED_HEADS $(collect(FUSE_EFIT_GATED_HEADS))")

    gstep = config["graphs"]["step"]
    gseq = config["graphs"]["sequence"]
    conv_shape = _shape_from_config(gstep["inputs"]["conv_state"], 1)
    ssm_shape = _shape_from_config(gstep["inputs"]["ssm_state"], 1)
    length(conv_shape) == 4 || error("fuse_efit: unexpected conv_state rank $(length(conv_shape)) (expected 4)")
    length(ssm_shape) == 5 || error("fuse_efit: unexpected ssm_state rank $(length(ssm_shape)) (expected 5)")
    seq_len = Int(gseq["seq_len"])
    period_s = Float64(config["cadence"]["period_ms"]) / 1000

    # the two gate thresholds: one sub-dict of outputs.hmode_thresholds carries `balanced` and `conservative`
    thr = nothing
    for v in values(config["outputs"]["hmode_thresholds"])
        if v isa AbstractDict && haskey(v, "balanced") && haskey(v, "conservative") && v["balanced"] isa Real
            thr = v
            break
        end
    end
    thr === nothing && error("fuse_efit: no balanced/conservative H-mode thresholds in $(bdir)/model_config.json")

    ranges_path = joinpath(bdir, "input_ranges.json")
    input_ranges = isfile(ranges_path) ? JSON.parsefile(ranges_path; dicttype=Dict{String,Any}) : nothing

    step = _load_inference_capped(joinpath(bdir, String(gstep["file"])))
    ORT = ONNXRunTime
    ORT.input_names(step) == collect(FUSE_EFIT_STEP_INPUTS) ||
        error("fuse_efit: step graph inputs $(ORT.input_names(step)) != $(collect(FUSE_EFIT_STEP_INPUTS))")
    ORT.output_names(step) == collect(FUSE_EFIT_STEP_OUTPUTS) ||
        error("fuse_efit: step graph outputs $(ORT.output_names(step)) != $(collect(FUSE_EFIT_STEP_OUTPUTS))")

    seq = Ref{Union{Nothing,ONNXRunTime.InferenceSession}}(nothing)
    return FuseEfitNN(step, seq, config, input_names, indicator_names, head_names, input_ranges, seq_len, form,
                      Float64(thr["balanced"]), Float64(thr["conservative"]), period_s, FUSE_EFIT_RHO_REF,
                      conv_shape, ssm_shape, String(bdir))
end

function _fuse_efit_seq_session(nn::FuseEfitNN)
    if nn.seq[] === nothing
        sess = _load_inference_capped(joinpath(nn.bundle_dir, String(nn.config["graphs"]["sequence"]["file"])))
        ONNXRunTime.input_names(sess) == collect(FUSE_EFIT_SEQ_INPUTS) ||
            error("fuse_efit: sequence graph inputs $(ONNXRunTime.input_names(sess)) != $(collect(FUSE_EFIT_SEQ_INPUTS))")
        nn.seq[] = sess
    end
    return nn.seq[]::ONNXRunTime.InferenceSession
end

# ── Inference ───────────────────────────────────────────────────────────────

"""
    fuse_efit_init_state(nn::FuseEfitNN) -> PedestalNNState

Zero state: the start of a shot (the exact initial condition the model was trained with).
"""
fuse_efit_init_state(nn::FuseEfitNN) = PedestalNNState(zeros(Float32, reverse(nn.conv_shape)), zeros(Float32, reverse(nn.ssm_shape)))

PedestalNNTracker(nn::FuseEfitNN, t_start::Float64) = PedestalNNTracker(fuse_efit_init_state(nn), FUSE_EFIT_N_HEADS, t_start)

"""
    fuse_efit_step(nn::FuseEfitNN, x, indicators, available, s::PedestalNNState) -> (pred::Vector{Float32}, s_new::PedestalNNState)

One 50 ms tick. `x` is the 25-vector of inputs in physical units
(`FUSE_EFIT_INPUT_NAMES` order), `indicators` the 5-vector of per-tick 0/1 flags
(`FUSE_EFIT_INDICATOR_NAMES` order) and `available` the 25-vector of per-shot
0/1 flags. Returns the 9 heads in physical units (`FUSE_EFIT_HEAD_NAMES` order)
and the advanced state. Functional: `s` is left untouched.

`x` may hold NaN in columns whose indicator or availability flag is 0 (they
cannot reach the model); NaN in a column that is in use is read by the graph as
the training mean, which [`build_fuse_efit_inputs`](@ref) uses on purpose for a
missing squareness value and nowhere else.
"""
function fuse_efit_step(nn::FuseEfitNN, x::AbstractVector{<:Real}, indicators::AbstractVector{<:Real},
                        available::AbstractVector{<:Real}, s::PedestalNNState)
    length(x) == FUSE_EFIT_N_INPUTS || error("fuse_efit_step: expected $(FUSE_EFIT_N_INPUTS) inputs, got $(length(x))")
    length(indicators) == FUSE_EFIT_N_INDICATORS || error("fuse_efit_step: expected $(FUSE_EFIT_N_INDICATORS) indicators, got $(length(indicators))")
    length(available) == FUSE_EFIT_N_INPUTS || error("fuse_efit_step: expected $(FUSE_EFIT_N_INPUTS) availability flags, got $(length(available))")
    sess = nn.step
    xb = Array{Float32}(undef, FUSE_EFIT_N_INPUTS, 1)          # reversed (1, 25)
    xb[:, 1] .= Float32.(x)
    ib = Array{Float32}(undef, FUSE_EFIT_N_INDICATORS, 1)      # reversed (1, 5)
    ib[:, 1] .= Float32.(indicators)
    ab = Array{Float32}(undef, FUSE_EFIT_N_INPUTS, 1)          # reversed (1, 25)
    ab[:, 1] .= Float32.(available)
    inputs = ONNXRunTime.CAPI.OrtValue[
        _ort_tensor(sess, xb, (1, FUSE_EFIT_N_INPUTS)),
        _ort_tensor(sess, ib, (1, FUSE_EFIT_N_INDICATORS)),
        _ort_tensor(sess, ab, (1, FUSE_EFIT_N_INPUTS)),
        _ort_tensor(sess, s.conv, nn.conv_shape),
        _ort_tensor(sess, s.ssm, nn.ssm_shape)]
    GC.@preserve xb ib ab s inputs begin
        outs = ONNXRunTime.CAPI.Run(sess.api, sess.session, nothing,
                                    collect(FUSE_EFIT_STEP_INPUTS), inputs, collect(FUSE_EFIT_STEP_OUTPUTS))
        pred = vec(_ort_output_native(sess, outs[1]))       # reversed (1, 9) -> 9-vector
        conv = _ort_output_native(sess, outs[2])
        ssm = _ort_output_native(sess, outs[3])
    end
    size(conv) == size(s.conv) && size(ssm) == size(s.ssm) || error("fuse_efit_step: unexpected state shapes $(size(conv)), $(size(ssm))")
    return pred, PedestalNNState(conv, ssm)
end

function _fuse_efit_check_shapes(fn::AbstractString, X::AbstractMatrix, I::AbstractMatrix, available::AbstractVector)
    size(X, 2) == FUSE_EFIT_N_INPUTS || error("$fn: expected inputs (K, $(FUSE_EFIT_N_INPUTS)), got $(size(X))")
    size(I) == (size(X, 1), FUSE_EFIT_N_INDICATORS) || error("$fn: expected indicators ($(size(X, 1)), $(FUSE_EFIT_N_INDICATORS)), got $(size(I))")
    length(available) == FUSE_EFIT_N_INPUTS || error("$fn: expected $(FUSE_EFIT_N_INPUTS) availability flags, got $(length(available))")
    return nothing
end

"""
    fuse_efit_rollout(nn::FuseEfitNN, X::AbstractMatrix, I::AbstractMatrix, available=ones(25)) -> Matrix{Float32}

Drive the step graph over a whole shot — inputs `(K, 25)`, indicators `(K, 5)` —
from the zero state. Returns `(K, 9)`. No length limit.
"""
function fuse_efit_rollout(nn::FuseEfitNN, X::AbstractMatrix{<:Real}, I::AbstractMatrix{<:Real},
                           available::AbstractVector{<:Real}=ones(Float32, FUSE_EFIT_N_INPUTS))
    _fuse_efit_check_shapes("fuse_efit_rollout", X, I, available)
    K = size(X, 1)
    out = Matrix{Float32}(undef, K, FUSE_EFIT_N_HEADS)
    s = fuse_efit_init_state(nn)
    for k in 1:K
        pred, s = fuse_efit_step(nn, @view(X[k, :]), @view(I[k, :]), available, s)
        out[k, :] .= pred
    end
    return out
end

"""
    fuse_efit_sequence(nn::FuseEfitNN, X::AbstractMatrix, I::AbstractMatrix, available=ones(25)) -> Matrix{Float32}

Score a whole `(K, 25)` shot with the fixed-length sequence graph. Shorter shots
are right-padded (zeros, indicators 1) and trimmed, which is exact because the
model is strictly causal; `K > nn.seq_len` is an error — use
[`fuse_efit_rollout`](@ref). The sequence graph is loaded on first use.
"""
function fuse_efit_sequence(nn::FuseEfitNN, X::AbstractMatrix{<:Real}, I::AbstractMatrix{<:Real},
                            available::AbstractVector{<:Real}=ones(Float32, FUSE_EFIT_N_INPUTS))
    _fuse_efit_check_shapes("fuse_efit_sequence", X, I, available)
    K = size(X, 1)
    K <= nn.seq_len || error("fuse_efit_sequence: shot has $K ticks but the sequence graph is fixed at $(nn.seq_len); use fuse_efit_rollout")
    xp = zeros(Float32, 1, nn.seq_len, FUSE_EFIT_N_INPUTS)
    xp[1, 1:K, :] .= Float32.(X)
    ip = ones(Float32, 1, nn.seq_len, FUSE_EFIT_N_INDICATORS)
    ip[1, 1:K, :] .= Float32.(I)
    ap = reshape(Float32.(collect(available)), 1, FUSE_EFIT_N_INPUTS)
    out = _fuse_efit_seq_session(nn)(Dict{String,Any}("inputs" => xp, "indicators" => ip, "available" => ap))
    pred = collect(Float32, out["predictions"])
    return pred[1, 1:K, :]
end

"""
    fuse_efit_named(pred::AbstractVector) -> NamedTuple

The 9-vector of heads as a NamedTuple keyed by `FUSE_EFIT_HEAD_NAMES`.
"""
fuse_efit_named(pred::AbstractVector) = NamedTuple{Symbol.(FUSE_EFIT_HEAD_NAMES)}(Tuple(Float64.(pred)))

"""
    fuse_efit_prediction(nn::FuseEfitNN, pred::AbstractVector, is_h_mode::Bool; time=NaN, n_ticks=0) -> NamedTuple

Assemble the `ActorPedestal.nn_prediction` record from one tick's raw heads and
the (hysteretic) H-mode gate decision. Carries the same fields the rest of the
actor reads from a fuse29 prediction (`ne`, `te_ped`, `ti_ped`, `t_rot_ped`,
`hmode_prob`, `is_h_mode`, `raw`, `rho_ref`, `time`, `n_ticks`) plus `model`.
The four H-mode-only heads (`FUSE_EFIT_GATED_HEADS`) are `NaN` unless
`is_h_mode`: they were trained on H-mode ticks only and are undefined in L-mode.
"""
function fuse_efit_prediction(nn::FuseEfitNN, pred::AbstractVector, is_h_mode::Bool; time::Float64=NaN, n_ticks::Int=0)
    p = Float64.(pred)
    gate(name, v) = (name in FUSE_EFIT_GATED_HEADS && !is_h_mode) ? NaN : v
    vals = Tuple(gate(name, p[i]) for (i, name) in enumerate(FUSE_EFIT_HEAD_NAMES))
    heads = NamedTuple{Symbol.(FUSE_EFIT_HEAD_NAMES)}(vals)
    return (; heads..., hmode_prob=p[end], is_h_mode, raw=Float32.(pred), rho_ref=nn.rho_ref, time, n_ticks, model=:fuse_efit)
end

# A column is "in use" when the shot has the channel and every indicator that governs it is 1.
function _fuse_efit_in_use(nn::FuseEfitNN, indicators::AbstractVector{<:Real}, available::AbstractVector{<:Real})
    used = [a > 0.5 for a in available]
    gov = nn.config["inputs"]["inputs"]["zeroed_when_indicator_is_0"]
    for (j, name) in enumerate(nn.input_names)
        for ind in gov[name]
            k = findfirst(==(String(ind)), nn.indicator_names)
            k === nothing && error("fuse_efit: unknown governing indicator `$ind` for channel `$name`")
            indicators[k] > 0.5 || (used[j] = false)
        end
    end
    return used
end

"""
    fuse_efit_check_units(nn::FuseEfitNN, x, indicators, available=ones(25); span_factor=nothing) -> Vector{String}

The model authors' units smoke test on one tick: a value in a column that is in
use and lies more than `span_factor` (3, from `input_ranges.json`) training
spans outside the training [p0.1, p99.9] range is reported. Catches amperes for
mega-amperes, centimetres for metres, watts for megawatts — errors that would
otherwise return confident numbers. Not a distribution-shift detector.
"""
function fuse_efit_check_units(nn::FuseEfitNN, x::AbstractVector{<:Real}, indicators::AbstractVector{<:Real},
                               available::AbstractVector{<:Real}=ones(Float32, FUSE_EFIT_N_INPUTS);
                               span_factor::Union{Nothing,Real}=nothing)
    nn.input_ranges === nothing && return ["input_ranges.json not present in this bundle; cannot check units"]
    factor = span_factor === nothing ? Float64(get(get(nn.input_ranges, "units_check_rule", Dict{String,Any}()), "span_factor", 3.0)) : Float64(span_factor)
    used = _fuse_efit_in_use(nn, indicators, available)
    ch = nn.input_ranges["channels"]
    warnings = String[]
    for (j, name) in enumerate(nn.input_names)
        (used[j] && isfinite(x[j])) || continue
        ref = ch[name]
        lo, hi = Float64(ref["pct"]["0.1"]), Float64(ref["pct"]["99.9"])
        width = hi - lo
        v = Float64(x[j])
        if v < lo - factor * width || v > hi + factor * width
            push!(warnings, "$name = $(round(v; sigdigits=4)) $(ref["unit"]) is more than $(factor) training spans outside the training p0.1..p99.9 range [$(round(lo; sigdigits=4)), $(round(hi; sigdigits=4))]. Check units.")
        end
    end
    return warnings
end

"""
    fuse_efit_in_distribution(nn::FuseEfitNN, x, indicators, available=ones(25)) -> Vector{NamedTuple}

Per-channel excursions of one input vector outside the training p1..p99
envelope, for the columns that are in use. A bookkeeping aid (the geometry of a
simulated equilibrium is the first place to look when fuse_efit and fuse29
disagree), not a validity test.
"""
function fuse_efit_in_distribution(nn::FuseEfitNN, x::AbstractVector{<:Real}, indicators::AbstractVector{<:Real},
                                   available::AbstractVector{<:Real}=ones(Float32, FUSE_EFIT_N_INPUTS); tol_frac::Float64=1e-3)
    nn.input_ranges === nothing && error("fuse_efit: input_ranges.json not present in $(nn.bundle_dir)")
    used = _fuse_efit_in_use(nn, indicators, available)
    ch = nn.input_ranges["channels"]
    excursions = NamedTuple{(:channel, :value, :lo, :hi),Tuple{String,Float64,Float64,Float64}}[]
    for (j, name) in enumerate(nn.input_names)
        (used[j] && isfinite(x[j])) || continue
        ref = ch[name]
        lo, hi = Float64(ref["pct"]["1"]), Float64(ref["pct"]["99"])
        pad = tol_frac * max(hi - lo, 1e-30)
        v = Float64(x[j])
        if v < lo - pad || v > hi + pad
            push!(excursions, (channel=name, value=v, lo=lo, hi=hi))
        end
    end
    return excursions
end
