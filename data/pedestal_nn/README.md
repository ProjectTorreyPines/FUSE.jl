# Pedestal NN data assets (fuse29)

Small assets for the fuse29 pedestal-predictor pipeline in
`src/actors/pedestal/`. The ONNX model weights deliberately **do not live
here** — they are distributed via HuggingFace and resolved at runtime via the
`FUSE_PEDESTAL_NN_DIR` environment variable (see
[`nn_predictor.jl`](../../src/actors/pedestal/nn_predictor.jl)).

## Contents

| File | Size | Purpose |
|---|---|---|
| `fuse29_reference_shot.npz` | ~24 KB | One real DIII-D discharge (shot 199055, test split) with the model's expected outputs for seed 0: `actuators` (K, 29) raw units, `predictions` (K, 9), `labels`/`label_mask`, `seed`, `period_ms`. Replayed by `test/test_pedestal_nn_predictor.jl` to verify routing, units, column order, timing and state handling of the Julia wrapper against the model authors' reference (tolerance 1e-3). Copied verbatim from the model repo's `examples/reference_shot.npz`. |

## The model

[`SCS-Lab/FUSE29-Pedestal-model`](https://huggingface.co/SCS-Lab/FUSE29-Pedestal-model)
(source: [NoahH72/FUSE-Pedestal-Predictor-Mamba-ONNX](https://github.com/NoahH72/FUSE-Pedestal-Predictor-Mamba-ONNX))
is a 6-layer Mamba-2 state-space model for DIII-D:

- **Inputs**: 29 actuator channels in **raw physical units** every 50 ms —
  `pinj` (MW), `tinj` (N·m), `ech_total` (MW), the 18 shaping F-coil currents
  `f1a..f9b` (A), `ecoila`/`ecoilb` (A), five gas valves `gasa_cal..gase_cal`
  (Torr·L/s), `bt` (T, signed). Normalisation is inside the graph. No shot
  history, no aux features, no plasma-response inputs (`ip`, `pohm` are gone).
- **Outputs**: `ne`, `te_ped`, `ti_ped`, `t_rot_ped` — profile values at
  ρ_tor = 0.85, valid in every regime; `neped_prmtan`, `teped_prmtan`,
  `rho_sym`, `ne_top_loc` — tanh-fit pedestal height/location, **H-mode only**;
  `hmode` — H-mode probability.
- **Stateful**: `fuse29_step.onnx` carries a ~3.1 MB recurrent state that must
  be zero at shot start and advanced exactly once per 50 ms of simulated time.
  `ActorPedestal` keeps it in `dd._aux[:fuse29]` (`Fuse29Tracker`) and handles
  catch-up, zero-order hold, iterated steps and rewinds.

HuggingFace layout (what `FUSE_PEDESTAL_NN_DIR` must contain):

```
manifest.json
s0/  fuse29_step.onnx  fuse29_seq.onnx  model_config.json  norms.json  actuator_ranges.json  provenance.json
s1/ … s7/   (optional: independently trained seeds, for spread across training runs)
```

Seed 0 is the champion and the default (`FUSE_PEDESTAL_NN_SEED` selects another).
One seed is ~90 MB; all eight are ~724 MB.

## Why ONNX weights are *not* in this repo

1. **Size** — two ~45 MB graphs per seed. Committing them would permanently
   bloat FUSE.jl's git history.
2. **Canonical home** — published on HuggingFace, SHA-addressed
   (`manifest.json` carries per-file sha256).
3. **Env-var contract** — `resolve_pedestal_nn_dir` prefers
   `ENV["FUSE_PEDESTAL_NN_DIR"]`, so containers bind-mount or bake in the
   artefact directory. Default when unset: the sibling checkout
   `<FUSE>/../pedestal-predictor-mamba-onnx-fuse/artifacts/`.

## Auto-fetch is on by default

A fresh `using FUSE; FUSE.load_pedestal_nn()` on a clean machine notices the
bundle is not there and pulls **seed 0** from HuggingFace on the spot — a
stdlib-only `Downloads.download` loop, no Python tooling.

```julia
using FUSE
nn = FUSE.load_pedestal_nn()          # ~90 MB on first call, idempotent on repeats
nn = FUSE.load_pedestal_nn(; seed=3)  # another seed (fetched on demand)
```

### Production containers — pre-bake at image build time

```dockerfile
ARG PEDESTAL_NN_HF_REVISION=main
ENV FUSE_PEDESTAL_NN_DIR=/opt/fuse29/artifacts
RUN julia --project=$FUSE_DIR -e \
    'using FUSE; FUSE.download_pedestal_nn!(ENV["FUSE_PEDESTAL_NN_DIR"]; \
                                            revision = ENV["PEDESTAL_NN_HF_REVISION"])'
```

After the build step the runtime auto-fetch hook is a no-op (it checks for
completeness before attempting any network call).

### Bind-mount mode (no auto-fetch)

```bash
export FUSE_PEDESTAL_NN_DIR=/nfs/fuse29/artifacts
export FUSE_PEDESTAL_NN_AUTOFETCH=0
```

Accepted opt-out values: `0`, `false`, `no`, `off` (case-insensitive).

### Pinning a revision

`download_pedestal_nn!` defaults to `revision="main"`. Production containers
should pin a HuggingFace commit SHA via the `revision` kwarg or
`ENV["FUSE_PEDESTAL_NN_HF_REVISION"]`. The repo id is overridable via
`ENV["FUSE_PEDESTAL_NN_HF_REPO"]` for in-house mirrors.

## Environment variables

| Variable | Default | Meaning |
|---|---|---|
| `FUSE_PEDESTAL_NN_DIR` | `<FUSE>/../pedestal-predictor-mamba-onnx-fuse/artifacts` | artefact root (or a bare bundle dir) |
| `FUSE_PEDESTAL_NN_SEED` | `0` | which `s<seed>/` bundle to load |
| `FUSE_PEDESTAL_NN_AUTOFETCH` | on | `0/false/no/off` disables the HuggingFace fetch |
| `FUSE_PEDESTAL_NN_HF_REPO` | `SCS-Lab/FUSE29-Pedestal-model` | HuggingFace repo id |
| `FUSE_PEDESTAL_NN_HF_REVISION` | `main` | HuggingFace revision / commit SHA |
| `FUSE_ONNX_THREADS` | `min(nthreads, 8)` | ONNX Runtime intra-op thread cap |

## Using it in a simulation

```julia
act.ActorPedestal.ne_from = :nn_predictor     # selects fuse29
act.ActorPedestal.fpe_source = :dd            # or :zmq (GSLite via ActorZMQ)
act.ActorPedestal.nn_ped_quantities = :ne_lh  # or :all (Te/Ti/rotation from the NN too)
act.ActorPedestal.nn_hmode_enter = 0.7        # Schmitt trigger on the hmode probability
act.ActorPedestal.nn_hmode_exit = 0.3
act.ActorPedestal.nn_warmup_time = 1.0        # optional: replay actuators from t=1 s before the first prediction
```

Things to know (from the model's own documentation):

- **Units have no guard rail.** `pinj` is MW, coils are A. The actor runs
  `fuse29_check_units` once on the first live actuator vector and warns.
- **Median fill.** Channels the host cannot supply (e.g. gas valves on the
  `:dd` path) are held at the training-split median, never the mean or zero.
- **Warm-up.** Outputs need ~1 s (20 ticks) from the zero state to settle.
- **Long steady state.** Under constant actuators the `hmode` probability decays
  slowly beyond the ~12.8 s training window; the hysteretic gate absorbs this,
  but treat the gate with suspicion in very long steady-state runs.
- **`ne`/`te_ped`/`ti_ped` are values at ρ_tor = 0.85**, not pedestal-top
  values at 0.9; the actor pins the profiles at 0.85.

## Tests

```bash
cd FUSE
julia --project test/test_pedestal_nn_predictor.jl   # contract, reference-shot replay, tracker semantics
julia --project test/smoke_pedestal_nn_shot.jl       # synthetic shot through the ZMQ actuator path
```
