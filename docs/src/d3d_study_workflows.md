# Getting Started with FUSE

FUSE (Fusion Unified Simulation Environment) is a Julia-based integrated plasma simulation framework. This guide walks you through setup, launching Jupyter, and running time-dependent DIII-D simulations using two complementary workflows:

- **Predictive** — forecasts plasma evolution from a given set of input parameters. Requires coupling with GSLite to provide those inputs, making it particularly suited for real-time control and digital twin applications.
- **Postdictive** — uses experimental data (kinetic profiles, equilibrium) from a real shot as input and evolves the plasma forward in time with physics models. This is the general-purpose workflow for validating FUSE against experiment.

---

## 1. Setup

You have two options:

### Option A — Local Installation (recommended if you want to modify the code)

Follow the official installation guide: https://fuse.help/dev/install.html

This installs Julia and FUSE on your laptop and gives you full access to the source code.

### Option B — Omega Server (no setup required)

Connect to the GA Omega cluster and use the pre-installed public version of FUSE.

Full instructions: https://fuse.help/dev/install_omega.html


---

## 2. Loading FUSE in a Notebook

Once Jupyter is open, create a new notebook and run:

```julia
using Plots
@time using FUSE
```

The first load may take a moment. After that you are ready to run simulations.

---

## 3. Postdictive Simulation

A postdictive simulation uses experimental data (kinetic profiles, equilibrium) from a real shot as input and evolves the plasma forward in time using physics models. It is the primary workflow for validating FUSE transport models against experiment, comparing simulated profiles with measured diagnostics across the full discharge. No external coupling is required — any DIII-D shot number can be used directly.

**A few things to note:**
- The first time you load a new shot, it will take a while to download the experimental data. With `use_local_cache=true` (default), subsequent runs will be much faster.
- `sty.reconstruction` controls which mode to use:
  - `true` — reconstruction mode: all profiles come directly from experiment, no transport solver
  - `false` — interpredictive mode: includes a transport solver (default: GKNN). If you have FUSE installed locally, you can switch to TGLFNN by setting `act.ActorTGLF.model = :TGLFNN` in `FUSE/src/studies/experiment_postdictive.jl`.

```julia
using Distributed

sty = FUSE.study_parameters(:Postdictive)
sty.server = "localhost"
sty.n_workers = 0
sty.release_workers_after_run = true
sty.save_folder = mktempdir("/your_path")
sty.device = :D3D
sty.shots = [200000]
sty.reconstruction = true           # true = reconstruction, false = interpredictive

sty.kw_case_parameters = Dict{Symbol,Any}(
    :use_local_cache => true,
    :fit_profiles => true,
    :EFIT_tree => "EFIT02")

@everywhere import FUSE
@everywhere ProgressMeter = FUSE.ProgressMeter

study = FUSE.StudyPostdictive(sty)
FUSE.run(study)
```

Results are saved to a timestamped subfolder inside `your_path/`. The full path is printed in the output, for example:
```
[ Info: saving simulation results to: /your_path/jl_5aKDmB/D3D_200000__2026-05-07T11:34:19.093__2230
```

You can load and inspect the results:

```julia
dd_sim = IMAS.json2imas("/your_path/.../dd_sim.json")
dd_exp = IMAS.json2imas("/your_path/.../dd_exp.json")

# Compare simulation vs experiment at t = 2.5 s
plot(dd_sim.core_profiles, time0=2.5)
plot!(dd_exp.core_profiles, time0=2.5)
```

---

## 4. Predictive Simulation

A predictive simulation starts from plasma parameters (device geometry, heating, targets) and predicts the plasma state using density predictor and transport models. This is the mode used for coupling with the digital twin through GSLite. The default transport solver in predictive mode is FINN.

Note that the full `dd` structure is not saved in predictive mode, in order to reduce runtime and storage. Instead, two output files are produced:
- `traces.json` — time traces of `li` and `beta_pol` from both FUSE and EFIT
- `validation_residuals.json` — validation residuals comparing time derivatives between FUSE and EFIT

These can be loaded the same way as shown above.

```julia
using Distributed

sty = FUSE.study_parameters(:PredictiveRT)
sty.server = "localhost"
sty.n_workers = 0
sty.release_workers_after_run = true
sty.save_folder = mktempdir("/your_path")
sty.device = :D3D
sty.shots = [200000]
sty.start_time = 0.5  # EFIT01 can be unreliable to use in FUSE early in the discharge; delay the start if needed
sty.end_time = 5.5

sty.kw_case_parameters = Dict{Symbol,Any}(
    :use_local_cache => true,
    :fit_profiles => true,
    :EFIT_tree => "EFIT01")

@everywhere import FUSE
@everywhere ProgressMeter = FUSE.ProgressMeter

study = FUSE.StudyPredictiveRT(sty)
@time FUSE.run(study)
```

---

## 5. Useful Links

- Documentation: https://fuse.help
- Installation: https://fuse.help/dev/install.html
- Omega setup: https://fuse.help/dev/install_omega.html
