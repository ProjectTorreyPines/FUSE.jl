# Pedestal NN data assets

Small binary assets consumed by the pedestal-predictor NN pipeline in
`src/actors/pedestal/`. Large ONNX model weights deliberately **do not live
here** — they are distributed via HuggingFace and resolved at runtime via the
`FUSE_PEDESTAL_NN_DIR` environment variable (see
[`nn_predictor.jl`](../../src/actors/pedestal/nn_predictor.jl)).

## Contents

| File | Size | Purpose |
|---|---|---|
| `history_morning.npz` | ~90 KB | DIII-D MSE history snapshot (50 shots × 446 features, shots 206717–206804) produced by `pedestal-predictor-onnx/scripts/dump_history_npz.py` on an Omega-cluster MDSplus run. Used as the default `FUSE_PEDESTAL_NN_HISTORY_NPZ` for smoke tests and real-time containers. |

The NPZ schema matches `HistoryManager.save_state` — see
[`nn_predictor.jl::load_history_npz`](../../src/actors/pedestal/nn_predictor.jl) for the
key-by-key contract.

## Why we ship a pre-computed NPZ (source-of-truth split)

The pedestal predictor has **two** distinct data inputs, fed from two
different places at runtime:

| Input | Shape | Source | Update cadence |
|---|---|---|---|
| 32-channel FPE actuator sequence | `(T, 32)` | **GSLite** via ZMQ → `dd._aux[:zmq_*]` → `build_fpe_sequences_from_aux` | every control cycle (real-time) |
| 50-slot × 446/458-feature MSE history | `(50, 446 or 458)` | **This NPZ** → `load_history_npz` → `dd._aux[:nn_history_buffer]` | operator-refreshed (not live) |

GSLite only provides real-time actuator/state signals; it does **not**
compute or supply the 50-shot machine-state history window the MSE encoder
needs. That window requires fetching ~67 PTDATA signals from MDSplus for
each of the last 50 shots and rolling them up into the feature vector —
an offline job that doesn't fit GSLite's realtime wire contract.

For the real-time container, the practical consequence is:

1. **Live path** — `ActorZMQ.receive!` keeps the `:zmq_*` aux records
   fresh from GSLite every cycle. Accurate per time step.
2. **History path** — a pre-computed NPZ (committed here) is loaded once at
   FUSE startup and held constant for the container's entire session.
   Updated only when an operator runs `dump_history_npz.py` against
   MDSplus and pushes a new NPZ into this directory.

This is acceptable today because the MSE history captures a *slow* machine
context (boronization state, last-week coil usage, typical gas mix) that
is stable across a run-day. A future enhancement (`compute_shot_stats` in
`nn_predictor.jl`, currently stubbed) would let FUSE compute this shot's
stats directly from `dd` and append to the buffer via
`push_shot_history!`, removing the MDSplus dependency entirely — but that
still doesn't live on the GSLite wire, it would be FUSE-side bookkeeping.

## Why ONNX weights are *not* in this repo

1. **Size** — 854 MB across 5 bundles, with `fpe_encoder.onnx` at 112 MB
   per bundle (over GitHub's 100 MB per-file hard limit). Committing them
   would permanently bloat FUSE.jl's git history.
2. **Canonical home** — models are already published on HuggingFace at
   [`SCS-Lab/pedestal-predictor-onnx`](https://huggingface.co/SCS-Lab/pedestal-predictor-onnx),
   SHA-addressed and cached.
3. **Env-var contract** — `resolve_pedestal_nn_dir` already prefers
   `ENV["FUSE_PEDESTAL_NN_DIR"]`, so containers just bind-mount or bake in
   a `onnx_models/` directory.

## Auto-fetch is on by default

A fresh `using FUSE; FUSE.load_pedestal_nn()` on a clean machine will
notice the bundles aren't there and pull them from
[`SCS-Lab/pedestal-predictor-onnx`](https://huggingface.co/SCS-Lab/pedestal-predictor-onnx)
on the spot — no Python tooling, no manual setup. The downloader is a
stdlib-only `Downloads.download` loop driven by the file manifest in
[`download_pedestal_nn!`](../../src/actors/pedestal/nn_predictor.jl).

This is intentional: the operator UX is "just run FUSE."

```julia
using FUSE
nn = FUSE.load_pedestal_nn()  # ~700 MB on first call, idempotent on repeats
```

### Production containers — pre-bake at image build time

Auto-fetch on first GSLite handshake means a one-time ~700 MB download
during the first prediction. Pre-bake at build time to keep first-cycle
latency tight:

```dockerfile
ARG PEDESTAL_NN_HF_REVISION=main
ENV FUSE_PEDESTAL_NN_DIR=/opt/pedestal-onnx/onnx_models
RUN julia --project=$FUSE_DIR -e \
    'using FUSE; FUSE.download_pedestal_nn!(ENV["FUSE_PEDESTAL_NN_DIR"]; \
                                            revision = ENV["PEDESTAL_NN_HF_REVISION"])'

# Runtime: point FUSE at the history NPZ shipped with this repo
ENV FUSE_PEDESTAL_NN_HISTORY_NPZ=/opt/FUSE/data/pedestal_nn/history_morning.npz
```

After the build step the directory is fully populated, so the runtime
auto-fetch hook becomes a no-op (it checks for completeness before
attempting any network call).

### Bind-mount mode (no auto-fetch)

If you mount the ONNX dir from an NFS share or a read-only volume, opt
out of auto-fetch so the loader fails loudly on missing files instead of
trying to write into the mount:

```bash
export FUSE_PEDESTAL_NN_DIR=/nfs/pedestal-onnx/onnx_models
export FUSE_PEDESTAL_NN_AUTOFETCH=0
```

Accepted opt-out values: `0`, `false`, `no`, `off` (case-insensitive).

### Pinning a revision

`download_pedestal_nn!` defaults to `revision="main"` for ergonomics.
Production containers should pin a specific HuggingFace commit SHA for
reproducible image rebuilds — override via the `revision` kwarg or
`ENV["FUSE_PEDESTAL_NN_HF_REVISION"]`. The repo URL itself is overridable
via `ENV["FUSE_PEDESTAL_NN_HF_REPO"]` for in-house mirrors.

## Refreshing the NPZ

The NPZ is intentionally static at runtime — GSLite does not resupply it.
When you want to advance the reference window (e.g. between run-days, or
to pin a specific campaign's machine state), regenerate it offline from
MDSplus on a network-reachable host. See
[`pedestal-predictor-onnx/docs/MAC_SETUP.md`](https://github.com/NoahH72/pedestal-predictor-onnx/blob/main/docs/MAC_SETUP.md)
for the full walkthrough.

```bash
conda activate pedestal-onnx
python scripts/dump_history_npz.py \
    --history-shots 206717:206804 \
    --out /path/to/FUSE/data/pedestal_nn/history_morning.npz
```

Then commit the replaced NPZ. The file is small enough that committing a
new revision occasionally is cheap, and SHA-pinning the NPZ alongside the
ONNX weights gives reproducible real-time container builds.
