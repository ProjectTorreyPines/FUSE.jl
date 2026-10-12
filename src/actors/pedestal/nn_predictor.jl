#= ============================================================================ =#
#  fuse29 — DIII-D Mamba-2 pedestal predictor (ONNX)                             #
#                                                                                #
#  Published at https://huggingface.co/SCS-Lab/FUSE29-Pedestal-model             #
#  Source repo: https://github.com/NoahH72/FUSE-Pedestal-Predictor-Mamba-ONNX    #
#                                                                                #
#  29 DIII-D actuator channels in RAW PHYSICAL UNITS in, 9 pedestal heads in     #
#  physical units out, one tick per 50 ms. Normalisation is baked into the       #
#  graphs; there is no shot-history, aux or machine-state input. The model is    #
#  strictly causal and stateful: `fuse29_step.onnx` carries a recurrent state    #
#  (conv_state (6,B,1536,3) + ssm_state (6,B,8,128,128), ~3.1 MB) that must be   #
#  zero at shot start and advanced exactly once per 50 ms of simulated time.     #
#  `fuse29_seq.onnx` evaluates the same weights over a whole shot at a fixed     #
#  256-tick length (offline scoring / tests).                                    #
#                                                                                #
#  Heads (in order):                                                             #
#    ne, te_ped, ti_ped, t_rot_ped        profile values at rho_tor = 0.85,      #
#                                         valid in every confinement regime      #
#    neped_prmtan, teped_prmtan,          tanh-fit pedestal height / location,   #
#    rho_sym, ne_top_loc                  ONLY DEFINED IN H-MODE (gate on hmode) #
#    hmode                                H-mode probability (sigmoid applied)   #
#                                                                                #
#  Runtime: ONNXRunTime.jl (Microsoft onnxruntime C library, opset 17).          #
#                                                                                #
#  Artefacts live *outside* FUSE (~90 MB per seed). Their root directory is      #
#  resolved from, in order:                                                      #
#    (1) the `onnx_dir` kwarg passed to `load_pedestal_nn`,                      #
#    (2) the environment variable `FUSE_PEDESTAL_NN_DIR`,                        #
#    (3) `<pedestal-predictor-mamba-onnx-fuse repo>/artifacts/` alongside FUSE.  #
#  Missing bundles are auto-fetched from HuggingFace (seed 0 only by default).   #
#= ============================================================================ =#

import ONNXRunTime
import JSON
import Downloads

const PEDESTAL_NN_ENV = "FUSE_PEDESTAL_NN_DIR"
const PEDESTAL_NN_AUTOFETCH_ENV = "FUSE_PEDESTAL_NN_AUTOFETCH"
const PEDESTAL_NN_HF_REPO_ENV = "FUSE_PEDESTAL_NN_HF_REPO"
const PEDESTAL_NN_HF_REVISION_ENV = "FUSE_PEDESTAL_NN_HF_REVISION"
const PEDESTAL_NN_SEED_ENV = "FUSE_PEDESTAL_NN_SEED"

const PEDESTAL_NN_HF_REPO_DEFAULT = "SCS-Lab/FUSE29-Pedestal-model"
const PEDESTAL_NN_HF_REVISION_DEFAULT = "main"
const PEDESTAL_NN_DEFAULT_SEED = 0

# ── The fuse29 I/O contract (mirrors fuse29/spec.py; asserted against the bundle at load) ──
const FUSE29_ACTUATOR_NAMES = (
    "pinj", "tinj", "ech_total",
    "f1a", "f2a", "f3a", "f4a", "f5a", "f6a", "f7a", "f8a", "f9a",
    "f1b", "f2b", "f3b", "f4b", "f5b", "f6b", "f7b", "f8b", "f9b",
    "ecoila", "ecoilb",
    "gasa_cal", "gasb_cal", "gasc_cal", "gasd_cal", "gase_cal",
    "bt")
const FUSE29_ACTUATOR_UNITS = (
    "MW", "N*m", "MW",
    "A", "A", "A", "A", "A", "A", "A", "A", "A",
    "A", "A", "A", "A", "A", "A", "A", "A", "A",
    "A", "A",
    "Torr*L/s", "Torr*L/s", "Torr*L/s", "Torr*L/s", "Torr*L/s",
    "T")
const FUSE29_HEAD_NAMES = ("ne", "te_ped", "ti_ped", "t_rot_ped",
                           "neped_prmtan", "teped_prmtan", "rho_sym", "ne_top_loc",
                           "hmode")
const FUSE29_HEAD_UNITS = ("1e19 m^-3", "keV", "keV", "krad/s",
                           "1e19 m^-3", "keV", "rho_tor", "rho_tor",
                           "probability")
const FUSE29_GATED_HEADS = ("neped_prmtan", "teped_prmtan", "rho_sym", "ne_top_loc")
const FUSE29_N_ACTUATORS = length(FUSE29_ACTUATOR_NAMES)
const FUSE29_N_HEADS = length(FUSE29_HEAD_NAMES)
const FUSE29_STEP_INPUTS = ("actuators", "conv_state", "ssm_state")
const FUSE29_STEP_OUTPUTS = ("predictions", "conv_state_out", "ssm_state_out")
const FUSE29_BUNDLE_FILES = ("fuse29_step.onnx", "fuse29_seq.onnx", "model_config.json",
                             "norms.json", "actuator_ranges.json", "provenance.json")
const _TOPLEVEL_FILES = ("manifest.json",)

"""
    Fuse29NN

One loaded fuse29 bundle (a single training seed): the stateful step graph, the
fixed-length sequence graph and the JSON sidecars.

Fields:
- `step`, `seq`        : ONNX Runtime sessions
- `config`             : parsed `model_config.json`
- `actuator_names`     : 29 input channel names, in graph column order
- `head_names`         : 9 output head names, in graph column order
- `actuator_ranges`    : parsed `actuator_ranges.json` (training-split distribution, or `nothing`)
- `seq_len`            : fixed sequence length of `fuse29_seq.onnx` (256)
- `seed`               : training seed of this bundle
- `period_s`           : model cadence (0.05 s)
- `rho_ref`            : radius the continuous heads are defined at (rho_tor = 0.85)
- `conv_shape`, `ssm_shape` : recurrent state shapes for batch 1
- `bundle_dir`         : directory the artefacts were loaded from
"""
struct Fuse29NN
    step::ONNXRunTime.InferenceSession
    seq::ONNXRunTime.InferenceSession
    config::Dict{String,Any}
    actuator_names::Vector{String}
    head_names::Vector{String}
    actuator_ranges::Union{Nothing,Dict{String,Any}}
    seq_len::Int
    seed::Int
    period_s::Float64
    rho_ref::Float64
    conv_shape::NTuple{4,Int}
    ssm_shape::NTuple{5,Int}
    bundle_dir::String
end

function Base.show(io::IO, nn::Fuse29NN)
    print(io, "Fuse29NN(seed=$(nn.seed), $(FUSE29_N_ACTUATORS)→$(FUSE29_N_HEADS), period=$(nn.period_s)s, rho_ref=$(nn.rho_ref), dir=\"$(nn.bundle_dir)\")")
end

"""
    Fuse29State

Recurrent state of the fuse29 step graph for a single shot (batch 1).
Zeros are the correct initial condition at shot start. `fuse29_step` never
mutates the state it is given; it returns a new one.

The arrays are stored in ONNX Runtime's native row-major layout, i.e. as Julia
column-major arrays with the ONNX dimensions **reversed** (`conv` is
`(3, 1536, 1, 6)` for the ONNX `(6, batch, 1536, 3)`, `ssm` is
`(128, 128, 8, 1, 6)` for `(6, batch, 8, 128, 128)`). `vec(arr)` is then exactly
the buffer ORT wants, so a tick costs no `permutedims` of the ~3 MB state — the
generic ONNXRunTime.jl path permutes it on the way in and out and was ~50x
slower than the graph itself.
"""
struct Fuse29State
    conv::Array{Float32,4}
    ssm::Array{Float32,5}
end
Base.copy(s::Fuse29State) = Fuse29State(copy(s.conv), copy(s.ssm))

# ── Path resolution ─────────────────────────────────────────────────────────

"""
    resolve_pedestal_nn_dir(; onnx_dir=nothing) -> String

Resolve the fuse29 artefact root: the directory that holds `manifest.json` and
the per-seed bundle subdirectories `s0/ … s7/` (a bare bundle directory that
itself contains `model_config.json` is also accepted by `load_pedestal_nn`).
"""
function resolve_pedestal_nn_dir(; onnx_dir::Union{Nothing,AbstractString}=nothing)
    if onnx_dir !== nothing
        return String(onnx_dir)
    end
    if haskey(ENV, PEDESTAL_NN_ENV)
        return ENV[PEDESTAL_NN_ENV]
    end
    return normpath(joinpath(__FUSE__, "..", "pedestal-predictor-mamba-onnx-fuse", "artifacts"))
end

function _default_seed()
    s = tryparse(Int, strip(get(ENV, PEDESTAL_NN_SEED_ENV, "")))
    return s === nothing ? PEDESTAL_NN_DEFAULT_SEED : s
end

# ── HuggingFace auto-fetch ──────────────────────────────────────────────────

"""
    pedestal_nn_files(; seed=0) -> Vector{String}

Relative paths (under the artefact root) that [`download_pedestal_nn!`](@ref)
fetches from HuggingFace for one seed: `manifest.json` plus the six files of
`s<seed>/`. Same set of files that `load_pedestal_nn` reads.
"""
function pedestal_nn_files(; seed::Integer=PEDESTAL_NN_DEFAULT_SEED)
    files = String[String(f) for f in _TOPLEVEL_FILES]
    for f in FUSE29_BUNDLE_FILES
        push!(files, "s$seed/$f")
    end
    return files
end

"""
    download_pedestal_nn!(target_dir::AbstractString;
                          repo=PEDESTAL_NN_HF_REPO_DEFAULT,
                          revision=PEDESTAL_NN_HF_REVISION_DEFAULT,
                          seed=0, force=false, verbose=true) -> target_dir

Pull one fuse29 seed bundle from HuggingFace into `target_dir`, producing the
directory layout `load_pedestal_nn` expects (`manifest.json` + `s<seed>/…`).
Idempotent: files already present on disk are skipped unless `force=true`.

`repo` / `revision` pin a specific HF mirror or commit SHA and are also
overridable via `ENV[$(repr(PEDESTAL_NN_HF_REPO_ENV))]` and
`ENV[$(repr(PEDESTAL_NN_HF_REVISION_ENV))]`. One seed is ~92 MB
(two ~45 MB ONNX graphs + JSON sidecars); all eight seeds are ~724 MB.

# Default behavior

`load_pedestal_nn` calls this automatically when the configured artefact root
is missing the requested seed. Opt out with
`ENV["$PEDESTAL_NN_AUTOFETCH_ENV"] = "0"` (or `"false"` / `"no"` / `"off"`)
when bind-mounting a read-only directory, so the loader fails fast instead of
trying to write into the mount.

# Container build usage (eager, pre-baked)

```
ARG PEDESTAL_NN_HF_REVISION=main
ENV FUSE_PEDESTAL_NN_DIR=/opt/fuse29/artifacts
RUN julia --project=\$FUSE_DIR -e \\
    'using FUSE; FUSE.download_pedestal_nn!(ENV["FUSE_PEDESTAL_NN_DIR"]; revision=ENV["PEDESTAL_NN_HF_REVISION"])'
```
"""
function download_pedestal_nn!(target_dir::AbstractString;
        repo::AbstractString=get(ENV, PEDESTAL_NN_HF_REPO_ENV, PEDESTAL_NN_HF_REPO_DEFAULT),
        revision::AbstractString=get(ENV, PEDESTAL_NN_HF_REVISION_ENV, PEDESTAL_NN_HF_REVISION_DEFAULT),
        seed::Integer=_default_seed(),
        force::Bool=false,
        verbose::Bool=true)
    mkpath(target_dir)
    files = pedestal_nn_files(; seed)
    base_url = "https://huggingface.co/$repo/resolve/$revision"
    verbose && @info "fuse29: fetching $(length(files)) files (seed $seed) from $repo @ $revision -> $target_dir"
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
            error("download_pedestal_nn!: failed to fetch $url\n  -> $err")
        end
        if filesize(dest) == 0
            rm(dest; force=true)
            error("download_pedestal_nn!: $rel downloaded empty from $url")
        end
        verbose && @info "  [$i/$(length(files))] -> $(filesize(dest)) B"
    end
    verbose && @info "fuse29: download complete ($(length(files)) files)"
    return String(target_dir)
end

"""
    pedestal_nn_dir_complete(dir::AbstractString; seed=0) -> Bool

True iff `dir` already contains every file [`pedestal_nn_files`](@ref) lists
for `seed` (zero-byte stragglers from interrupted downloads count as missing),
or `dir` is itself a complete bare bundle directory.
"""
function pedestal_nn_dir_complete(dir::AbstractString; seed::Integer=_default_seed())
    isdir(dir) || return false
    _bundle_complete(dir) && return true
    for rel in pedestal_nn_files(; seed)
        path = joinpath(dir, rel)
        (isfile(path) && filesize(path) > 0) || return false
    end
    return true
end

function _bundle_complete(bdir::AbstractString)
    isdir(bdir) || return false
    for f in FUSE29_BUNDLE_FILES
        path = joinpath(bdir, f)
        (isfile(path) && filesize(path) > 0) || return false
    end
    return true
end

# Auto-fetch is on by default so a fresh `using FUSE; load_pedestal_nn()`
# Just Works on a clean machine without manual setup.
function _maybe_autofetch!(dir::AbstractString; seed::Integer=_default_seed())
    autofetch = lowercase(strip(get(ENV, PEDESTAL_NN_AUTOFETCH_ENV, "")))
    autofetch in ("0", "false", "no", "off") && return false
    pedestal_nn_dir_complete(dir; seed) && return false
    @info "fuse29: ONNX bundle for seed $seed missing or incomplete; auto-fetching from HuggingFace ($(get(ENV, PEDESTAL_NN_HF_REPO_ENV, PEDESTAL_NN_HF_REPO_DEFAULT)) @ $(get(ENV, PEDESTAL_NN_HF_REVISION_ENV, PEDESTAL_NN_HF_REVISION_DEFAULT))) -> $dir\n  Set ENV[\"$PEDESTAL_NN_AUTOFETCH_ENV\"] = \"0\" to opt out."
    download_pedestal_nn!(dir; seed)
    return true
end

# ── Loading ─────────────────────────────────────────────────────────────────

"""
    load_pedestal_nn(; onnx_dir=nothing, seed=nothing) -> Fuse29NN

Load one fuse29 bundle. `onnx_dir` (or `ENV["FUSE_PEDESTAL_NN_DIR"]`) may be
either the artefact root holding `manifest.json` and `s0/ … s7/`, in which case
`seed` (default `ENV["FUSE_PEDESTAL_NN_SEED"]`, else 0 — the champion) selects
the bundle, or a bare bundle directory containing `model_config.json`.
Cached per process; repeated calls return the same object.

The bundle is validated against the contract in this file at load time
(actuator and head names/order, graph input/output names, state shapes), so a
mismatched or mislabelled bundle fails here rather than producing plausible
nonsense downstream.
"""
function load_pedestal_nn(; onnx_dir::Union{Nothing,AbstractString}=nothing,
                            seed::Union{Nothing,Integer}=nothing)
    dir = resolve_pedestal_nn_dir(; onnx_dir)
    s = seed === nothing ? _default_seed() : Int(seed)
    return _cached_pedestal_nn(dir, s)
end

const _PEDESTAL_NN_CACHE = Dict{Tuple{String,Int},Fuse29NN}()

function _cached_pedestal_nn(dir::AbstractString, seed::Int)
    cache_key = (abspath(dir), seed)
    haskey(_PEDESTAL_NN_CACHE, cache_key) && return _PEDESTAL_NN_CACHE[cache_key]
    root = cache_key[1]
    _bundle_complete(root) || _maybe_autofetch!(root; seed)
    isdir(root) || error("fuse29: artefact directory does not exist: $root\n" *
        "Auto-fetch was skipped (ENV[\"$PEDESTAL_NN_AUTOFETCH_ENV\"] is set to opt-out) " *
        "or failed before creating the directory. Either:\n" *
        "  • Unset / set ENV[\"$PEDESTAL_NN_AUTOFETCH_ENV\"]=\"1\" and retry (downloads from HuggingFace)\n" *
        "  • Run `FUSE.download_pedestal_nn!(\"$root\")` explicitly\n" *
        "  • Set ENV[\"$PEDESTAL_NN_ENV\"] or pass `onnx_dir=...` to point at an existing copy")
    bdir = _resolve_bundle_dir(root, seed)
    nn = _load_fuse29(bdir, seed)
    _PEDESTAL_NN_CACHE[cache_key] = nn
    return nn
end

# Mirrors fuse29.runtime.from_pretrained: a bare bundle dir wins, else `s<seed>/` under the root.
function _resolve_bundle_dir(root::AbstractString, seed::Int)
    isfile(joinpath(root, "model_config.json")) && return String(root)
    bdir = joinpath(root, "s$seed")
    _bundle_complete(bdir) && return bdir
    missing_files = [f for f in FUSE29_BUNDLE_FILES if !(isfile(joinpath(bdir, f)) && filesize(joinpath(bdir, f)) > 0)]
    error("fuse29: bundle directory $bdir is missing $(missing_files). " *
          "Expected a per-seed bundle (`$(root)/s$seed/`) or a bare bundle directory with model_config.json.")
end

# ---------------------------------------------------------------------------
# ONNX Runtime thread cap.
#
# `ONNXRunTime.load_inference` creates the session with default options, so the
# intra-op pool gets one thread per *host* core with per-core affinity. Inside a
# container / Slurm cpuset that is fatal: on a 128-core omega node with a 32-CPU
# cpuset every pthread_setaffinity_np outside the set fails (EINVAL, ~1000 log
# lines) and the session hung the coupled run; on a 256-thread Perlmutter node
# it is a 256-thread pool per model for a handful of tiny inferences. Cap the
# pool (FUSE_ONNX_THREADS, default min(Threads.nthreads(), 8)); with an explicit
# thread count ORT does not apply affinities.
# ---------------------------------------------------------------------------
const ONNX_THREADS_ENV = "FUSE_ONNX_THREADS"

function _ort_threads()
    n = tryparse(Int, get(ENV, ONNX_THREADS_ENV, ""))
    return n === nothing ? max(1, min(Threads.nthreads(), 8)) : max(1, n)
end

# separate function: an interpolated-pointer @ccall inside try/catch trips Julia's lowering
function _set_intra_op_threads!(api, session_options, n::Integer)
    status = @ccall $(api.SetIntraOpNumThreads)(session_options::Ptr{Cvoid}, Cint(n)::Cint)::Ptr{Cvoid}
    ONNXRunTime.CAPI.check_and_release(api, status)
    return nothing
end

function _load_inference_capped(path::AbstractString; threads::Int=_ort_threads())
    ORT = ONNXRunTime
    try
        api = ORT.GetApi(; execution_provider=:cpu)
        env = ORT.CAPI.CreateEnv(api; name="defaultenv", logging_level=ORT.parse_logging_level(:warning))
        so = ORT.CAPI.CreateSessionOptions(api)
        _set_intra_op_threads!(api, so, threads)
        session = ORT.CAPI.CreateSession(api, env, path, so)
        meminfo = ORT.CAPI.CreateCpuMemoryInfo(api)
        allocator = ORT.CAPI.CreateAllocator(api, session, meminfo)
        ins = ORT.input_names(api, session, allocator)
        outs = ORT.output_names(api, session, allocator)
        @debug "fuse29: ONNX session for $(basename(path)) with $threads intra-op thread(s)"
        return ORT.InferenceSession(api, :cpu, session, meminfo, allocator, ins, outs)
    catch e
        @warn "fuse29: could not create a thread-capped ONNX session (ONNXRunTime.jl API changed?); falling back to load_inference with default threading" exception = e
        return ORT.load_inference(path)
    end
end

function _shape_from_config(dims::AbstractVector, batch::Int)
    return Tuple(Int[d == "batch" ? batch : Int(d) for d in dims])
end

function _load_fuse29(bdir::AbstractString, seed::Int)
    config = JSON.parsefile(joinpath(bdir, "model_config.json"); dicttype=Dict{String,Any})

    # -- contract assertions (names & order are load-bearing) --
    actuator_names = String.(config["inputs"]["names"])
    head_names = String.(config["outputs"]["names"])
    actuator_names == collect(FUSE29_ACTUATOR_NAMES) ||
        error("fuse29: actuator names/order in $(bdir)/model_config.json do not match FUSE29_ACTUATOR_NAMES:\n  bundle: $actuator_names\n  FUSE:   $(collect(FUSE29_ACTUATOR_NAMES))")
    head_names == collect(FUSE29_HEAD_NAMES) ||
        error("fuse29: head names/order in $(bdir)/model_config.json do not match FUSE29_HEAD_NAMES:\n  bundle: $head_names\n  FUSE:   $(collect(FUSE29_HEAD_NAMES))")
    gated = haskey(config["outputs"], "hmode_gated") ? String.(config["outputs"]["hmode_gated"]["heads"]) : collect(FUSE29_GATED_HEADS)
    gated == collect(FUSE29_GATED_HEADS) ||
        error("fuse29: hmode-gated heads in bundle ($gated) differ from FUSE29_GATED_HEADS $(collect(FUSE29_GATED_HEADS))")

    gstep = config["graphs"]["fuse29_step.onnx"]
    gseq = config["graphs"]["fuse29_seq.onnx"]
    conv_shape = _shape_from_config(gstep["inputs"]["conv_state"], 1)
    ssm_shape = _shape_from_config(gstep["inputs"]["ssm_state"], 1)
    length(conv_shape) == 4 || error("fuse29: unexpected conv_state rank $(length(conv_shape)) (expected 4)")
    length(ssm_shape) == 5 || error("fuse29: unexpected ssm_state rank $(length(ssm_shape)) (expected 5)")
    seq_len = Int(get(gseq, "seq_len", gseq["inputs"]["actuators"][2]))

    cadence = get(config, "cadence", Dict{String,Any}())
    period_s = Float64(get(cadence, "period_ms", 50.0)) / 1000
    rho_ref = Float64(get(cadence, "rho_tor", 0.85))
    bundle_seed = Int(get(config, "seed", seed))
    bundle_seed == seed || @warn "fuse29: bundle in $bdir reports seed $bundle_seed but was requested as seed $seed"

    ranges_path = joinpath(bdir, "actuator_ranges.json")
    actuator_ranges = isfile(ranges_path) ? JSON.parsefile(ranges_path; dicttype=Dict{String,Any}) : nothing

    step = _load_inference_capped(joinpath(bdir, "fuse29_step.onnx"))
    seq = _load_inference_capped(joinpath(bdir, "fuse29_seq.onnx"))
    ORT = ONNXRunTime
    ORT.input_names(step) == collect(FUSE29_STEP_INPUTS) ||
        error("fuse29: fuse29_step.onnx inputs $(ORT.input_names(step)) != $(collect(FUSE29_STEP_INPUTS))")
    ORT.output_names(step) == collect(FUSE29_STEP_OUTPUTS) ||
        error("fuse29: fuse29_step.onnx outputs $(ORT.output_names(step)) != $(collect(FUSE29_STEP_OUTPUTS))")
    ORT.input_names(seq) == ["actuators"] ||
        error("fuse29: fuse29_seq.onnx inputs $(ORT.input_names(seq)) != [\"actuators\"]")

    return Fuse29NN(step, seq, config, actuator_names, head_names, actuator_ranges,
                    seq_len, bundle_seed, period_s, rho_ref, conv_shape, ssm_shape, String(bdir))
end

# ── Inference ───────────────────────────────────────────────────────────────

"""
    fuse29_init_state(nn::Fuse29NN) -> Fuse29State

Zero state: the start of a shot. Zeros are the correct initial condition, not a
placeholder — they are what the sequence graph implicitly assumes before tick 0.
"""
fuse29_init_state(nn::Fuse29NN) = Fuse29State(zeros(Float32, reverse(nn.conv_shape)), zeros(Float32, reverse(nn.ssm_shape)))

# Wrap a Julia array that already holds ORT's row-major buffer (reversed dims) as an OrtValue of `onnx_shape`.
function _ort_tensor(sess::ONNXRunTime.InferenceSession, data::Array{Float32}, onnx_shape)
    size(data) == reverse(Tuple(onnx_shape)) || error("fuse29: array of size $(size(data)) is not the reversed ONNX shape $(Tuple(onnx_shape))")
    return ONNXRunTime.CAPI.CreateTensorWithDataAsOrtValue(sess.api, sess.meminfo, vec(data), onnx_shape)
end

# Copy an output OrtValue into an owned Julia array with reversed ONNX dims (no permutation).
function _ort_output_native(sess::ONNXRunTime.InferenceSession, val)
    GC.@preserve val begin
        pda = ONNXRunTime.CAPI.unsafe_GetTensorMutableData(sess.api, val)
        return copy(parent(pda))
    end
end

"""
    fuse29_step(nn::Fuse29NN, u::AbstractVector{<:Real}, s::Fuse29State) -> (pred::Vector{Float32}, s_new::Fuse29State)

One 50 ms tick. `u` is the 29-vector of actuators in raw physical units
(`FUSE29_ACTUATOR_NAMES` order). Returns the 9 heads in physical units
(`FUSE29_HEAD_NAMES` order) and the advanced state. Functional: `s` is left
untouched, so a host that iterates inside a time step can re-step from the same
committed state any number of times.

Goes through the ONNX Runtime C API directly with the state in ORT's native
layout (see [`Fuse29State`](@ref)); ~5 ms per tick on one core.
"""
function fuse29_step(nn::Fuse29NN, u::AbstractVector{<:Real}, s::Fuse29State)
    length(u) == FUSE29_N_ACTUATORS || error("fuse29_step: expected $(FUSE29_N_ACTUATORS) actuators, got $(length(u))")
    all(isfinite, u) || error("fuse29_step: non-finite actuator value(s): $(collect(zip(FUSE29_ACTUATOR_NAMES, u))[.!isfinite.(u)])")
    sess = nn.step
    x = Array{Float32}(undef, FUSE29_N_ACTUATORS, 1)      # reversed (1, 29)
    x[:, 1] .= Float32.(u)
    inputs = ONNXRunTime.CAPI.OrtValue[
        _ort_tensor(sess, x, (1, FUSE29_N_ACTUATORS)),
        _ort_tensor(sess, s.conv, nn.conv_shape),
        _ort_tensor(sess, s.ssm, nn.ssm_shape)]
    GC.@preserve x s inputs begin
        outs = ONNXRunTime.CAPI.Run(sess.api, sess.session, nothing,
                                    collect(FUSE29_STEP_INPUTS), inputs, collect(FUSE29_STEP_OUTPUTS))
        pred = vec(_ort_output_native(sess, outs[1]))       # reversed (1, 9) -> 9-vector
        conv = _ort_output_native(sess, outs[2])
        ssm = _ort_output_native(sess, outs[3])
    end
    size(conv) == size(s.conv) && size(ssm) == size(s.ssm) || error("fuse29_step: unexpected state shapes $(size(conv)), $(size(ssm))")
    return pred, Fuse29State(conv, ssm)
end

"""
    fuse29_rollout(nn::Fuse29NN, U::AbstractMatrix{<:Real}) -> Matrix{Float32}

Drive the step graph over a whole `(K, 29)` actuator sequence from the zero
state. Returns `(K, 9)`. No length limit.
"""
function fuse29_rollout(nn::Fuse29NN, U::AbstractMatrix{<:Real})
    size(U, 2) == FUSE29_N_ACTUATORS || error("fuse29_rollout: expected (K, $(FUSE29_N_ACTUATORS)), got $(size(U))")
    K = size(U, 1)
    out = Matrix{Float32}(undef, K, FUSE29_N_HEADS)
    s = fuse29_init_state(nn)
    for k in 1:K
        pred, s = fuse29_step(nn, @view(U[k, :]), s)
        out[k, :] .= pred
    end
    return out
end

"""
    fuse29_sequence(nn::Fuse29NN, U::AbstractMatrix{<:Real}) -> Matrix{Float32}

Score a whole `(K, 29)` shot with the fixed-length sequence graph. Shorter shots
are right-padded with zeros and trimmed (exact, since the model is strictly
causal); `K > nn.seq_len` is an error — use [`fuse29_rollout`](@ref).
"""
function fuse29_sequence(nn::Fuse29NN, U::AbstractMatrix{<:Real})
    size(U, 2) == FUSE29_N_ACTUATORS || error("fuse29_sequence: expected (K, $(FUSE29_N_ACTUATORS)), got $(size(U))")
    K = size(U, 1)
    K <= nn.seq_len || error("fuse29_sequence: shot has $K ticks but the sequence graph is fixed at $(nn.seq_len); use fuse29_rollout")
    padded = zeros(Float32, 1, nn.seq_len, FUSE29_N_ACTUATORS)
    padded[1, 1:K, :] .= Float32.(U)
    out = nn.seq(Dict{String,Any}("actuators" => padded))
    pred = collect(Float32, out["predictions"])
    return pred[1, 1:K, :]
end

"""
    fuse29_named(pred::AbstractVector) -> NamedTuple

The 9-vector of heads as a NamedTuple keyed by `FUSE29_HEAD_NAMES`.
"""
fuse29_named(pred::AbstractVector) = NamedTuple{Symbol.(FUSE29_HEAD_NAMES)}(Tuple(Float64.(pred)))

"""
    fuse29_prediction(nn::Fuse29NN, pred::AbstractVector, is_h_mode::Bool; time=NaN, n_ticks=0) -> NamedTuple

Assemble the `ActorPedestal.nn_prediction` record from one tick's raw heads and
the (hysteretic) H-mode gate decision. The four H-mode-only heads
(`FUSE29_GATED_HEADS`) are `NaN` unless `is_h_mode`, on purpose: they were never
trained on L-mode ticks and the graph emits unconstrained numbers for them there.
"""
function fuse29_prediction(nn::Fuse29NN, pred::AbstractVector, is_h_mode::Bool; time::Float64=NaN, n_ticks::Int=0)
    p = Float64.(pred)
    gate(name, v) = (name in FUSE29_GATED_HEADS && !is_h_mode) ? NaN : v
    vals = Tuple(gate(name, p[i]) for (i, name) in enumerate(FUSE29_HEAD_NAMES))
    heads = NamedTuple{Symbol.(FUSE29_HEAD_NAMES)}(vals)
    return (; heads..., hmode_prob=p[end], is_h_mode, raw=Float32.(pred), rho_ref=nn.rho_ref, time, n_ticks)
end

"""
    fuse29_median_actuators(nn::Fuse29NN) -> Vector{Float32}

Training-split median of every channel (from `actuator_ranges.json`), the
recommended stand-in for a channel the host cannot supply. Not the mean: the
mean vector is jointly unphysical (`bt` is signed and bimodal, ECH is zero on
most ticks), and zeros put the coil channels several sigma out.
"""
function fuse29_median_actuators(nn::Fuse29NN)
    nn.actuator_ranges === nothing && error("fuse29: actuator_ranges.json not present in $(nn.bundle_dir); cannot build the median actuator vector")
    ch = nn.actuator_ranges["channels"]
    return Float32[Float64(ch[name]["pct"]["50"]) for name in nn.actuator_names]
end

"""
    fuse29_check_units(nn::Fuse29NN, U; tol_sigma=12.0) -> Vector{String}

Compare actuators (a 29-vector or a `(K, 29)` matrix, raw units) against the
training distribution. Catches the one preprocessing mistake that cannot raise:
watts where the graph wants megawatts, kA where it wants A. Returns
human-readable warnings; empty means every channel is within `tol_sigma` of the
training mean. A smoke test for unit errors, not a distribution-shift detector.
"""
function fuse29_check_units(nn::Fuse29NN, U::AbstractVecOrMat{<:Real}; tol_sigma::Float64=12.0)
    nn.actuator_ranges === nothing && return ["actuator_ranges.json not present in this bundle; cannot check units"]
    X = U isa AbstractVector ? reshape(Float64.(U), 1, :) : Float64.(U)
    size(X, 2) == FUSE29_N_ACTUATORS || error("fuse29_check_units: expected 29 columns, got $(size(X, 2))")
    warnings = String[]
    ch = nn.actuator_ranges["channels"]
    for (j, name) in enumerate(nn.actuator_names)
        ref = ch[name]
        std = max(Float64(ref["std"]), 1e-12)
        col = @view X[:, j]
        z = maximum(abs.((col .- Float64(ref["mean"])) ./ std))
        if z > tol_sigma
            push!(warnings, "$name: observed [$(round(minimum(col); sigdigits=4)), $(round(maximum(col); sigdigits=4))] reaches $(round(Int, z)) sigma from the training mean $(round(Float64(ref["mean"]); sigdigits=4)). Training p1..p99 was [$(round(Float64(ref["pct"]["1"]); sigdigits=4)), $(round(Float64(ref["pct"]["99"]); sigdigits=4))] $(ref["units"]). Check units.")
        end
    end
    return warnings
end

"""
    fuse29_in_distribution(nn::Fuse29NN, u::AbstractVector; tol_frac=1e-3) -> Vector{NamedTuple}

Per-channel excursions of one actuator vector outside the training p1..p99
envelope (each bound padded by `tol_frac` of the channel's own span). A
bookkeeping aid for annotating counterfactual scans, not a validity test — the
check is per channel and cannot see whether the *combination* has ever been run.
"""
function fuse29_in_distribution(nn::Fuse29NN, u::AbstractVector{<:Real}; tol_frac::Float64=1e-3)
    nn.actuator_ranges === nothing && error("fuse29: actuator_ranges.json not present in $(nn.bundle_dir)")
    ch = nn.actuator_ranges["channels"]
    excursions = NamedTuple{(:channel, :value, :lo, :hi),Tuple{String,Float64,Float64,Float64}}[]
    for (j, name) in enumerate(nn.actuator_names)
        ref = ch[name]
        lo, hi = Float64(ref["pct"]["1"]), Float64(ref["pct"]["99"])
        pad = tol_frac * max(hi - lo, Float64(ref["std"]), 1e-30)
        v = Float64(u[j])
        if v < lo - pad || v > hi + pad
            push!(excursions, (channel=name, value=v, lo=lo, hi=hi))
        end
    end
    return excursions
end

# ── Time-step bookkeeping for a host simulation ─────────────────────────────

"""
    Fuse29Tracker

Per-shot bookkeeping that guarantees the recurrent state is advanced exactly
once per 50 ms of simulated time, however often the host evaluates the pedestal
inside a step (implicit coupling, rejected/repeated steps). Lives in
`dd._aux[:fuse29]` so it survives actor re-creation.

- `committed*` : state/time/gate after the last *accepted* host step
- `pending*`   : state/time/gate advanced to the most recent `t_now` (not yet committed)
- `last`       : most recent 9-vector, for zero-order hold when the host steps faster than 50 ms
- `t_start`    : start of the model's 50 ms grid (shot start, or `nn_warmup_time`)
"""
mutable struct Fuse29Tracker
    committed::Fuse29State
    committed_time::Float64
    committed_hmode::Bool
    pending::Fuse29State
    pending_time::Float64
    pending_hmode::Bool
    t_start::Float64
    last::Vector{Float32}
    n_total::Int
    units_checked::Bool
end

function Fuse29Tracker(nn::Fuse29NN, t_start::Float64)
    s = fuse29_init_state(nn)
    return Fuse29Tracker(s, t_start, false, s, t_start, false, t_start, fill(NaN32, FUSE29_N_HEADS), 0, false)
end

"""
    fuse29_advance!(tr::Fuse29Tracker, nn::Fuse29NN, t_now::Float64, actuators_at::Function;
                    hmode_enter=0.7, hmode_exit=0.3) -> (pred::Vector{Float32}, is_h_mode::Bool, n_ticks::Int)

Bring the model up to `t_now`. `actuators_at(t)` must return the 29-vector of
raw-unit actuators for the 50 ms interval ending at `t`.

1. If `t_now` is past the pending time, the pending state is committed (exactly
   one commit per accepted host step). A re-run at the same `t_now` commits
   nothing and re-steps from the committed state, so iterating inside a step
   never advances the model's clock.
2. `n = floor((t_now - committed_time) / period)` ticks are stepped, each fed the
   actuators for its own 50 ms sub-interval (catch-up when the host steps slower
   than 50 ms; warm-up replay from `t_start`).
3. `n == 0` (host steps faster than 50 ms): zero-order hold of the last output.

The H-mode gate is a Schmitt trigger: enter H-mode when `hmode > hmode_enter`,
leave when `hmode < hmode_exit`, ignore the ambiguous middle — `hmode` is right
~94% of the time with the errors concentrated at L–H / H–L transitions, so a
bare threshold chatters there. Note also that under a constant actuator request
the `hmode` probability decays slowly beyond the ~12.8 s the model was trained
on (the state encodes elapsed shot time); the hysteresis makes the gate robust
to that, but very long steady-state runs should treat the gate with suspicion.
"""
function fuse29_advance!(tr::Fuse29Tracker, nn::Fuse29NN, t_now::Float64, actuators_at::Function;
                         hmode_enter::Float64=0.7, hmode_exit::Float64=0.3)
    period = nn.period_s
    eps = 1e-3 * period
    hmode_exit <= hmode_enter || error("fuse29_advance!: hmode_exit ($hmode_exit) must be <= hmode_enter ($hmode_enter)")
    t_now >= tr.committed_time - eps || error("fuse29_advance!: t_now=$t_now is before the committed time $(tr.committed_time); reset the tracker for a rewind")

    # commit exactly once per accepted host step
    if t_now > tr.pending_time + eps
        tr.committed = tr.pending
        tr.committed_time = tr.pending_time
        tr.committed_hmode = tr.pending_hmode
    end

    n = floor(Int, (t_now - tr.committed_time) / period + 1e-6)
    s = tr.committed
    h = tr.committed_hmode
    pred = tr.last
    for i in 1:n
        t_i = tr.committed_time + i * period
        u = actuators_at(t_i)
        pred, s = fuse29_step(nn, u, s)
        p = Float64(pred[end])
        h = h ? (p > hmode_exit) : (p > hmode_enter)
    end
    if n > 0
        tr.pending = s
        tr.pending_time = tr.committed_time + n * period
        tr.pending_hmode = h
        tr.last = pred
        tr.n_total += n
    end
    return pred, h, n
end
