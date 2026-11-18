# Pedestal NN data assets (fuse29 and fuse_efit)

Small assets for the two pedestal-predictor pipelines in
`src/actors/pedestal/`: **fuse29** (the default) and **fuse_efit** (its
successor with geometry inputs; see [the fuse_efit section](#fuse_efit) below). The ONNX model weights deliberately **do not live
here** — they are distributed via HuggingFace and resolved at runtime via the
`FUSE_PEDESTAL_NN_DIR` environment variable (see
[`nn_predictor.jl`](../../src/actors/pedestal/nn_predictor.jl)).

## Contents

| File | Size | Purpose |
|---|---|---|
| `fuse29_reference_shot.npz` | ~24 KB | One real DIII-D discharge (shot 199055, test split) with the model's expected outputs for seed 0: `actuators` (K, 29) raw units, `predictions` (K, 9), `labels`/`label_mask`, `seed`, `period_ms`. Replayed by `test/test_pedestal_nn_predictor.jl` to verify routing, units, column order, timing and state handling of the Julia wrapper against the model authors' reference (tolerance 1e-3). Copied verbatim from the model repo's `examples/reference_shot.npz`. |
| `fuse_efit_reference_shot.npz` | ~23 KB | One real DIII-D discharge (shot 193190, in-era validation split, never trained on) with the fuse_efit outputs for both forms: `inputs` (K, 25) physical units, `indicators` (K, 5), `available` (25), `expected_single` / `expected_averaged` (K, 9), `tolerance` (9). Replayed by `test/test_pedestal_nn_efit.jl`. Copied verbatim from `examples/reference_shot.npz` on the model repo's `fuse-efit` branch. |

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
act.ActorPedestal.ne_from = :nn_predictor     # the edge comes from a pedestal predictor
act.ActorPedestal.nn_model = :fuse29          # the default; :fuse_efit selects the other model (below)
act.ActorPedestal.fpe_source = :dd            # or :zmq (GSLite via ActorZMQ)
act.ActorPedestal.nn_ped_quantities = :ne_lh  # or :all (Te/Ti/rotation from the NN too)
act.ActorPedestal.nn_hmode_enter = 0.7        # Schmitt trigger on the hmode probability (these two are
act.ActorPedestal.nn_hmode_exit = 0.3         # the fuse29 defaults when left unset)
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
julia --project test/test_pedestal_nn_dd_actuators.jl  # :dd actuators, coil convention, injected torque
julia --project test/test_pedestal_nn_efit.jl          # fuse_efit: contract, reference-shot replay, masking, tracker
julia --project test/test_pedestal_nn_efit_inputs.jl   # fuse_efit: boundary-shape inputs on analytic equilibria, input builder
```

## fuse_efit

[`SCS-Lab/FUSE-PED`](https://huggingface.co/SCS-Lab/FUSE-PED) (upstream name
FUSE-EFIT; code and documentation on the `fuse-efit` branch of
[NoahH72/FUSE-Pedestal-Predictor-Mamba-ONNX](https://github.com/NoahH72/FUSE-Pedestal-Predictor-Mamba-ONNX):
`HANDOFF.md`, `docs/EXPORT.md`, `docs/INPUT_CONTRACT.md`, `MODEL_CARD.md`) is
the same kind of model — Mamba-2, one tick per 50 ms, the same nine heads at
ρ_tor = 0.85 — with a different input set. Both models stay available so that
validation runs can be made with each and compared.

| | fuse29 | fuse_efit |
|---|---|---|
| inputs | 29: `pinj`, `tinj`, `ech_total`, **20 coil currents**, 5 gas valves, `bt` | 25: `pinj`, `tinj`, `ech_mw`, 5 gas valves, `bt`, **`ip` [MA]**, **14 boundary-shape quantities**, **`vloop` [V]** |
| extra input tensors | – | `indicators` (5 flags per tick), `available` (25 flags per shot) |
| a channel the host cannot supply | held at the training median | flagged unavailable for the shot (a trained state) |
| recurrent state | ~3.1 MB | ~2.5 MB (`:single`), ~20 MB (`:averaged`) |
| H-mode-only heads | `neped_prmtan`, `teped_prmtan`, `rho_sym`, `ne_top_loc` | `ne_knee`, `te_knee`, `rho_sym`, `ne_top_loc` |
| default gate (enter / leave) | 0.7 / 0.3 | the form's validation thresholds: 0.54 / 0.50 (`:single`), 0.53 / 0.47 (`:averaged`) |
| forms | 8 independent seeds | `:single` (one network) and `:averaged` (eight averaged every tick) |

The 14 shape quantities are `aminor`, `rsurf`, `zsurf`, `kappa`, `tritop`,
`tribot`, `drsep` (clipped to ±0.1 m), `rvsout`, `gapin`, `rvsou`, `rvsiu`,
`shafranov`, `zeta_ou`, `zeta_ol`; the five indicators are `efit_valid`,
`vloop_valid` and one presence flag per strike point. FUSE computes them with
`FUSE.fuse_efit_shape` following the definitions of the model's
`docs/INPUT_CONTRACT.md` (Shafranov shift = magnetic-axis R − geometric-centre
R; Luce squareness of the two outer quadrants from the discrete boundary
vertices, the same construction as `IMAS.luce_squareness`).

**The model was trained on EFIT01**, the between-shot magnetics-only
reconstruction. Feeding it another equilibrium costs accuracy: the authors
measured a 31 % higher loss with real-time EFIT, mostly through the Shafranov
shift, the squareness and the strike points. That is the reason to keep fuse29
next to it, and what `nn_state_source` is for.

### Using it

```julia
act.ActorPedestal.ne_from = :nn_predictor
act.ActorPedestal.nn_model = :fuse_efit
act.ActorPedestal.nn_efit_form = :single        # or :averaged (~10x the cost per tick, 1.5 GB download)
act.ActorPedestal.fpe_source = :dd              # actuators: :dd or :zmq, as for fuse29
act.ActorPedestal.nn_state_source = :simulation # Ip, shape and loop voltage from FUSE's own equilibrium
                                                # (FRESCO, or GSLite's flux map when coupled through ZMQ);
                                                # :replay takes them from the experimental equilibrium instead
```

| input | `nn_state_source = :simulation` | `:replay` |
|---|---|---|
| `ip` | `dd._aux[:zmq_Ip_avg]` with `fpe_source = :zmq`, else `IMAS.Ip(cp1d, eqt)` from the current-diffusion solution (the equilibrium's `global_quantities.ip` when there are no current profiles) | `global_quantities.ip` of the experimental slices in the tick window |
| 14 shape channels | `fuse_efit_shape` of the simulation's equilibrium slice | mean of `fuse_efit_shape` over the experimental slices in (t − 50 ms, t]; `efit_valid = 0` when there is none (no hold-last) |
| `vloop` | `IMAS.vloop(cp1d, eqt)` from the current-diffusion solution | dψ_boundary/dt of the experimental equilibrium (the dd carries no measured loop voltage) |

Warm-up ticks before the simulation starts (`nn_warmup_time`) always read the
stored experimental traces and the experimental equilibrium. For `:replay` runs
that should match training as closely as possible, load the shot with
`EFIT_tree="EFIT01"` (the D3D case defaults to EFIT02).

### How close FUSE's inputs are to the training inputs

Shot 193190 is the model's reference shot: `fuse_efit_reference_shot.npz` holds the
inputs the model's own pipeline built for it. Loading the same shot in FUSE with
`EFIT_tree="EFIT01"` and building the inputs with `nn_state_source = :replay`
gives, over its 130 ticks (difference in units of the channel's training
standard deviation):

| channels | median difference | note |
|---|---|---|
| `aminor`, `rsurf`, `zsurf`, `kappa`, `tritop`, `tribot`, `gapin`, `shafranov` | < 0.001 σ | the X-point is taken as the boundary's lowest point, as EFIT does |
| `ip`, `drsep`, all five indicators | < 0.01 σ | strike-point flags agree on 96 % of ticks or more |
| `bt` | 0.03 σ | 2.5 % smaller than PTDATA `bt` (the equilibrium's vacuum field) |
| `rvsout` | 0.07 σ | 9 mm outboard of EFIT's |
| `vloop` | 0.08 σ | dψ_boundary/dt against the measured loop voltage, correlation 0.83 |
| `pinj` | 7 % high | the launched power in the data FUSE loads against PTDATA `bmspinj` (the training spread of `pinj` is set by outliers, so no σ is given) |
| `tinj` | 0.12 σ | beam-geometry estimate, 10 % high: the torque per watt is within 3 %, the rest is the power |
| `zeta_ou`, `zeta_ol` | 0.28 σ, 0.15 σ | the resolution bias of the discrete-vertex squareness (0.016): FUSE's traced boundary has 150-180 points, EFIT01's about 90 |
| `gasa_cal` … `gase_cal` | not available | the DIII-D data export has no gas-valve flows |

**Beam torque.** `tinj` is the *injected* torque (PTDATA `bmstinj`). With
`fpe_source = :dd` FUSE computes it from the beam geometry and the launched
power, Σ 2 P R_tang / v over the beams and their three energy components,
averaged over the 50 ms tick; the torque the beam actor deposits in
`core_sources` is used only when the `dd` has no beam data. Against PTDATA the
estimate is 5 % high on 193190 and 5 % low on 199055 over a run, which is the
difference in launched power on those shots (+6 % and −7 %). The beam actor's
torque is the deposited one, not the injected one, and was 12 % high on 199055
(16 % after 2.5 s). With the estimate in its place the rotation predicted after
2.5 s on 199055 went from 17.6 to 15.4 krad/s for fuse29 (measured 15.3) and
from 31.4 to 27.8 krad/s for fuse_efit on the replayed state.

**Upper strike points of a lower-single-null plasma.** When the upper X-point is
the secondary one, FUSE reports its two strike points whenever its separatrix is
within 0.4 m of the boundary at the midplane. EFIT01 does not. Over 89 DIII-D
shots (8,182 slices with a secondary upper X-point) it reports them on nearly
all slices when the separatrices are within 1 cm of each other, on about 80 % at
1 cm to 3 cm, on 60 % at 3 cm to 5 cm and almost never beyond 6 cm, and no
geometric quantity tried (separation, X-point height, distance to the wall)
predicts which slices on better than 81 %: EFIT follows the second separatrix
from the midplane, and whether its contour tracer gets past the X-point decides
the outcome. EFIT01 also loses the upper X-point itself on some slices when it
is within a few centimetres of the wall, and then stores a saturated `drsep`.
Where EFIT01 does find it, FUSE's X-point and `drsep` agree with it to 1 mm.
On DIII-D 199055 (`drsep` between −0.04 m and −0.07 m) FUSE had
`rvsou_present = 1` on 72 of 74 ticks and EFIT01 on 13; with EFIT01's flags and
`drsep` the rotation predicted after 2.5 s is 7 krad/s lower (20.7 krad/s,
measured 15.3).
The lower strike points of an upper-single-null plasma are not affected: EFIT01
reports them on every slice.

**Gas valves.** With `fpe_source = :dd` the experimental `dd` carries no gas
flows (the OMAS DIII-D mapping covers the valve hardware only), so fuse_efit is
told the five gas channels are unavailable and fuse29 holds them at the training
median. With `fpe_source = :zmq` they arrive from GSLite.

**Simulation geometry.** With `nn_state_source = :simulation` the shape is that
of FUSE's own equilibrium. When the equilibrium solver does not converge on a
step, FUSE keeps the previous equilibrium and the shape inputs lag by that
step; the plasma current does not (it comes from the current-diffusion
solution). On DIII-D 199055 this held the shape for 13 consecutive ticks of the
current ramp (0.44 s to 1.03 s).

**Simulated or measured state.** With `:simulation` the Shafranov shift and the
loop voltage are those of the simulated plasma, so they carry the simulation's
own pressure and resistivity, and the predicted pedestal feeds back on its own
inputs. On DIII-D 193190 the simulated flat-top Shafranov shift was 0.07 m
where EFIT01 has 0.06 m (0.036 m after 3.5 s), and the simulated loop voltage
0.22 V to 0.37 V where the replayed one is 0.34 V to 0.43 V. That was enough to
change the prediction: the `:simulation` run had an RMSE of 0.172 keV on
`te_ped` and 0.247 keV on `ti_ped`, the `:replay` run 0.070 keV and 0.123 keV.
Feeding the model the recorded `:simulation` inputs offline with only these two
channels replaced by the replayed values gave the same predictions as the full
set of `:replay` inputs; the
other shape channels, `ip` and `tinj` made no difference. Use
`nn_state_source = :replay` to judge the model apart from the simulation.

### Artefacts

```
manifest.json  SHA256SUMS
era/single/    fuse_efit_single_step.onnx    fuse_efit_single_sequence.onnx    model_config.json  norms.json  input_ranges.json  provenance.json
era/averaged/  fuse_efit_averaged_step.onnx  fuse_efit_averaged_sequence.onnx  …  (fetched only when `nn_efit_form = :averaged`)
```

| Variable | Default | Meaning |
|---|---|---|
| `FUSE_PEDESTAL_EFIT_NN_DIR` | `<FUSE>/../pedestal-predictor-mamba-onnx-fuse/artifacts/fuse_efit` | download root (or its `era/` folder, or a bare bundle dir) |
| `FUSE_PEDESTAL_EFIT_NN_HF_REPO` | `SCS-Lab/FUSE-PED` | HuggingFace repo id |
| `FUSE_PEDESTAL_EFIT_NN_HF_REVISION` | `main` | HuggingFace revision / commit SHA |
| `FUSE_PEDESTAL_NN_AUTOFETCH` | on | shared with fuse29: `0/false/no/off` disables the fetch |

`FUSE.load_pedestal_efit_nn(; form=:single)` fetches the requested form on first
use (~184 MB for `:single`) and checks every file's size against
`manifest.json`; `FUSE.download_pedestal_efit_nn!(dir; form)` does it eagerly.

### Comparing the two models on a shot

`StudyPredictiveRT` takes `nn_model`, `nn_state_source` and `nn_efit_form`. Each
run saves `pedestal_nn.json`: every model tick with the inputs the model was
given, its nine raw heads, the gate decision, and the experimental `ne`, `Te`,
`Ti` and rotation at ρ = 0.85 on the same time base, with the RMSE of each.

```julia
sty = FUSE.study_parameters(:PredictiveRT)
sty.device = :D3D; sty.shots = [199055]; sty.reconstruction = false; sty.save_folder = "runs"
sty.nn_model = :fuse29;    FUSE.run(FUSE.StudyPredictiveRT(sty))
sty.nn_model = :fuse_efit; FUSE.run(FUSE.StudyPredictiveRT(sty))
FUSE.compare_pedestal_nn_runs("runs/D3D_199055__…__fuse29", "runs/D3D_199055__…__fuse_efit")
FUSE.plot_pedestal_nn_comparison("runs/D3D_199055__…__fuse29", "runs/D3D_199055__…__fuse_efit"; savefig="D3D_199055_pedestal_nn.png")
```

`compare_pedestal_nn_runs` prints the RMSE of each head against the experiment
on the ticks both runs share; `plot_pedestal_nn_comparison` draws the same heads
against the experimental points, with each model's H-mode probability and gate.
For two runs of the same model, name them with `labels`, e.g.
`labels=("simulation", "replay")` for two `nn_state_source` settings.

The model tag is appended to each run's folder name. The same record is
available in memory as `FUSE.pedestal_nn_history(dd)` after any run.
